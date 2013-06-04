///////////////////////////////////////////////////////////////////////////////
// 
//  Copyright (2013) Alexander Stukowski
//
//  This file is part of OVITO (Open Visualization Tool).
//
//  OVITO is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  OVITO is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//  The matrix decomposition code has been taken from the book
//  Graphics Gems IV - Ken Shoemake, Polar AffineTransformation Decomposition. 
//
///////////////////////////////////////////////////////////////////////////////

#include <base/Base.h>
#include <base/linalg/AffineDecomposition.h>
#include <base/linalg/Matrix4.h>

namespace Ovito {

void decomp_affine(Matrix4& A, AffineDecomposition* parts);

/******************************************************************************
* Constructor that decomposes the matrix into its affine parts.
******************************************************************************/
AffineDecomposition::AffineDecomposition(const AffineTransformation& tm)
{
    // Extends matrix to 4x4.
	Matrix4 A(tm);
    decomp_affine(A, this);

	OVITO_ASSERT_MSG(std::abs(scaling.Q.dot(scaling.Q) - 1.0) <= FLOATTYPE_EPSILON, "AffineDecomposition", "Resulting quaternion is not normalized.");

	// The following code checks whether the decomposed parts give the original affine matrix.
#if defined(_DEBUG) && 0

	// Decompose 4x4 affine matrix A as TFRUK(U transpose), where t contains the
	// translation components, q contains the rotation R, u contains U, k contains
	// scale factors, and f contains the sign of the determinant.

	AffineTransformation S = AffineTransformation::scaling(scaling);
	AffineTransformation R = AffineTransformation::rotation(rotation);
	AffineTransformation T = AffineTransformation::translation(translation);
	AffineTransformation F(IDENTITY);	
	F(0,0) *= f; F(1,1) *= f; F(2,2) *= f;
	AffineTransformation product = T * F * R * S;
    for(size_t i=0; i<4; i++) {
    	if(!product[i].Equals(tm[i], FLOATTYPE_EPSILON)) {
    		qDebug() << "Original matrix: " << tm;
    		qDebug() << "Product of affine parts: " << product;
			OVITO_ASSERT_MSG(false, "AffineDecomposition(const AffineTransformation& tm)", "Could not decompose matrix to affine transformations.");
    	}
    }

#endif
}

enum QuatPart {X, Y, Z, W};

/** Fill out 3x3 matrix to 4x4 **/
#define mat_pad(A) (A(W,X)=A(X,W)=A(W,Y)=A(Y,W)=A(W,Z)=A(Z,W)=0,A(W,W)=1)

/** Copy nxn matrix A to C using "gets" for assignment **/
#define mat_copy(C,gets,A,n) {int i,j; for(i=0;i<n;i++) for(j=0;j<n;j++) C(i,j) gets (A(i,j));}

/** Copy transpose of nxn matrix A to C using "gets" for assignment **/
#define mat_tpose(AT,gets,A,n) {int i,j; for(i=0;i<n;i++) for(j=0;j<n;j++) AT(i,j) gets (A(j,i));}

/** Assign nxn matrix C the element-wise combination of A and B using "op" **/
#define mat_binop(C,gets,A,op,B,n) {int i,j; for(i=0;i<n;i++) for(j=0;j<n;j++) C(i,j) gets (A(i,j)) op (B(i,j));}

/** Multiply the upper left 3x3 parts of A and B to get AB **/
inline void mat_mult(const Matrix4& A, const Matrix4& B, Matrix4& AB)
{
    for(int i=0; i<3; i++) 
		for(int j=0; j<3; j++)
			AB(i,j) = A(i,0)*B(0,j) + A(i,1)*B(1,j) + A(i,2)*B(2,j);
}

/** Return dot product of length 3 vectors va and vb **/
inline FloatType vdot(const Vector4& va, const Vector4& vb)
{
    return (va[0]*vb[0] + va[1]*vb[1] + va[2]*vb[2]);
}

/** Return dot product of length 3 vectors va and vb **/
inline FloatType vdot(const Vector3& va, const Vector4& vb)
{
    return (va[0]*vb[0] + va[1]*vb[1] + va[2]*vb[2]);
}

/** Set v to cross product of length 3 vectors va and vb **/
inline void vcross(const Vector4& va, const Vector4& vb, Vector4& v)
{
    v[0] = va[1]*vb[2] - va[2]*vb[1];
    v[1] = va[2]*vb[0] - va[0]*vb[2];
    v[2] = va[0]*vb[1] - va[1]*vb[0];
}

/** Set v to cross product of length 3 vectors va and vb **/
inline void vcross(const Vector4& va, const Vector4& vb, Vector3& v)
{
    v[0] = va[1]*vb[2] - va[2]*vb[1];
    v[1] = va[2]*vb[0] - va[0]*vb[2];
    v[2] = va[0]*vb[1] - va[1]*vb[0];
}

/** Set MadjT to transpose of inverse of M times determinant of M **/
void adjoint_transpose(Matrix4& M, Matrix4& MadjT)
{
	Vector4 v = Vector4::Zero();
    vcross(M.row(1), M.row(2), v);
	MadjT.setRow(0, v);
    vcross(M.row(2), M.row(0), v);
	MadjT.setRow(1, v);
    vcross(M.row(0), M.row(1), v);
	MadjT.setRow(2, v);
}


/******* Quaternion Preliminaries *******/

/* Return conjugate of quaternion. */
inline Quaternion Qt_Conj(const Quaternion& q)
{
	return q.inverse();
}

/* Construct a unit quaternion from rotation matrix.  Assumes matrix is
 * used to multiply column vector on the left: vnew = mat vold.	 Works
 * correctly for right-handed coordinate system and right-handed rotations.
 * Translation and perspective components ignored. */
Quaternion Qt_FromMatrix(const Matrix4& mat)
{
    /* This algorithm avoids near-zero divides by looking for a large component
     * - first w, then x, y, or z.  When the trace is greater than zero,
     * |w| is greater than 1/2, which is as small as a largest component can be.
     * Otherwise, the largest diagonal entry corresponds to the largest of |x|,
     * |y|, or |z|, one of which must be larger than |w|, and at least 1/2. */
    Quaternion qu;
    FloatType tr, s;

    tr = mat(X,X) + mat(Y,Y)+ mat(Z,Z);
    if(tr >= 0.0) {
	    s = sqrt(tr + mat(W,W));
	    qu.w() = s*0.5;
	    s = 0.5 / s;
	    qu.x() = (mat(Z,Y) - mat(Y,Z)) * s;
	    qu.y() = (mat(X,Z) - mat(Z,X)) * s;
	    qu.z() = (mat(Y,X) - mat(X,Y)) * s;
	} else {
	    int h = X;
	    if(mat(Y,Y) > mat(X,X)) h = Y;
	    if(mat(Z,Z) > mat(h,h)) h = Z;
	    switch(h) {
#define caseMacro(i,j,k,I,J,K)  case I:	s = sqrt( (mat(I,I) - (mat(J,J)+mat(K,K))) + mat(W,W) ); qu.i() = s*0.5;	s = 0.5 / s;	qu.j() = (mat(I,J) + mat(J,I)) * s;		qu.k() = (mat(K,I) + mat(I,K)) * s;		qu.w() = (mat(K,J) - mat(J,K)) * s;		break
	    caseMacro(x,y,z,X,Y,Z);
	    caseMacro(y,z,x,Y,Z,X);
	    caseMacro(z,x,y,Z,X,Y);
	    }
	}
    if (mat(W,W) != 1.0) qu /= sqrt(mat(W,W));
    return qu;
}

/******* Decomp Auxiliaries *******/

/** Compute either the 1 or infinity norm of M, depending on tpose **/
inline FloatType mat_norm(const Matrix4& M, bool tpose)
{
    int i;
    FloatType sum, max;
    max = 0.0;
    for (i=0; i<3; i++) {
		if (tpose) sum = abs(M(0,i))+abs(M(1,i))+abs(M(2,i));
		else	   sum = abs(M(i,0))+abs(M(i,1))+abs(M(i,2));
		if (max<sum) max = sum;
    }
    return max;
}

FloatType norm_inf(const Matrix4& M) {return mat_norm(M, false);}
FloatType norm_one(const Matrix4& M) {return mat_norm(M, true);}

/** Return index of column of M containing maximum abs entry, or -1 if M=0 **/
int find_max_col(const Matrix4& M)
{
    FloatType abs, max;
    int i, j, col;
    max = 0.0; col = -1;
	for(i=0; i<3; i++) {
		for(j=0; j<3; j++) {
			abs = M(i,j); if (abs<0.0) abs = -abs;
			if (abs>max) {max = abs; col = j;}
		}
    }
    return col;
}

/** Setup u for Household reflection to zero all v components but first **/
void make_reflector(const Vector3& v, Vector3& u)
{
    float s = v.length();
    u[0] = v[0]; u[1] = v[1];
    u[2] = v[2] + ((v[2]<0.0) ? -s : s);
    s = sqrt(2.0/u.squaredLength());
    u[0] = u[0]*s; u[1] = u[1]*s; u[2] = u[2]*s;
}

/** Apply Householder reflection represented by u to column vectors of M **/
void reflect_cols(Matrix4& M, const Vector3& u)
{
    int i, j;
	for (i=0; i<3; i++) {
		FloatType s = u[0]*M(0,i) + u[1]*M(1,i) + u[2]*M(2,i);
		for (j=0; j<3; j++) M(j,i) -= u[j]*s;
    }
}
/** Apply Householder reflection represented by u to row vectors of M **/
void reflect_rows(Matrix4& M, const Vector3& u)
{
    int i, j;
    for (i=0; i<3; i++) {
		FloatType s = vdot(u, M.column(i));
		for (j=0; j<3; j++) M(i,j) -= u[j]*s;
    }
}

/** Find orthogonal factor Q of rank 1 (or less) M **/
void do_rank1(Matrix4& M, Matrix4& Q)
{
    Vector3 v1, v2;
    FloatType s;
    int col;
	Q = Matrix4::Identity();
    /* If rank(M) is 1, we should find a non-zero column in M */
    col = find_max_col(M);
    if (col<0) return; /* Rank is 0 */
    v1[0] = M(0,col); v1[1] = M(1,col); v1[2] = M(2,col);
    make_reflector(v1, v1); reflect_cols(M, v1);
    v2[0] = M(2,0); v2[1] = M(2,1); v2[2] = M(2,2);
    make_reflector(v2, v2); reflect_rows(M, v2);
    s = M(2,2);
    if (s<0.0) Q(2,2) = -1.0;
    reflect_cols(Q, v1); reflect_rows(Q, v2);
}

/** Find orthogonal factor Q of rank 2 (or less) M using adjoint transpose **/
void do_rank2(Matrix4& M, Matrix4& MadjT, Matrix4& Q)
{
    Vector3 v1, v2;
    FloatType w, x, y, z, c, s, d;
    int col;
    /* If rank(M) is 2, we should find a non-zero column in MadjT */
    col = find_max_col(MadjT);
    if (col<0) {do_rank1(M, Q); return;} /* Rank<2 */
    v1[0] = MadjT(0,col); v1[1] = MadjT(1,col); v1[2] = MadjT(2,col);
    make_reflector(v1, v1); reflect_cols(M, v1);
    vcross(M.row(0), M.row(1), v2);
    make_reflector(v2, v2); reflect_rows(M, v2);
    w = M(0,0); x = M(0,1); y = M(1,0); z = M(1,1);
    if (w*z>x*y) {
	c = z+w; s = y-x; d = sqrt(c*c+s*s); c = c/d; s = s/d;
	Q(0,0) = Q(1,1) = c; Q(0,1) = -(Q(1,0) = s);
    } else {
	c = z-w; s = y+x; d = sqrt(c*c+s*s); c = c/d; s = s/d;
	Q(0,0) = -(Q(1,1) = c); Q(0,1) = Q(1,0) = s;
    }
    Q(0,2) = Q(2,0) = Q(1,2) = Q(2,1) = 0.0; Q(2,2) = 1.0;
    reflect_cols(Q, v1); reflect_rows(Q, v2);
}

/* Polar Decomposition of 3x3 matrix in 4x4,
 * M = QS.  See Nicholas Higham and Robert S. Schreiber,
 * Fast Polar Decomposition of An Arbitrary Matrix,
 * Technical Report 88-942, October 1988,
 * Department of Computer Science, Cornell University.
 */
FloatType polar_decomp(Matrix4& M, Matrix4& Q, Matrix4& S)
{
    Matrix4 Mk, MadjTk, Ek;
    FloatType det, M_one, M_inf, MadjT_one, MadjT_inf, E_one, gamma, g1, g2;
    int i, j;

    mat_tpose(Mk,=,M,3);
    M_one = norm_one(Mk);  M_inf = norm_inf(Mk);
    do {
		adjoint_transpose(Mk, MadjTk);
		det = vdot(Mk.row(0), MadjTk.row(0));
		if(det==0.0) {
			do_rank2(Mk, MadjTk, Mk); 
			break;
		}
		MadjT_one = norm_one(MadjTk); MadjT_inf = norm_inf(MadjTk);
		gamma = sqrt(sqrt((MadjT_one*MadjT_inf)/(M_one*M_inf))/abs(det));
		g1 = gamma*0.5;
		g2 = 0.5/(gamma*det);
		mat_copy(Ek,=,Mk,3);
		mat_binop(Mk,=,g1*Mk,+,g2*MadjTk,3);
		mat_copy(Ek,-=,Mk,3);
		E_one = norm_one(Ek);
		M_one = norm_one(Mk);  M_inf = norm_inf(Mk);
    } while (E_one>(M_one*FLOATTYPE_EPSILON));
    mat_tpose(Q,=,Mk,3); 
	mat_pad(Q);
    mat_mult(Mk, M, S);	 
	mat_pad(S);
    for(i=0; i<3; i++) 
		for (j=i; j<3; j++)
			S(i,j) = S(j,i) = 0.5*(S(i,j)+S(j,i));
    return det;
}

/******* Spectral Decomposition *******/

/* Compute the spectral decomposition of symmetric positive semi-definite S.
 * Returns rotation in U and scale factors in result, so that if K is a diagonal
 * matrix of the scale factors, then S = U K (U transpose). Uses Jacobi method.
 * See Gene H. Golub and Charles F. Van Loan. AffineTransformation Computations. Hopkins 1983.
 */
Vector3 spect_decomp(Matrix4& S, Matrix4& U)
{
	Vector3 kv;
	FloatType Diag[3],OffD[3]; // OffD is off-diag (by omitted index)
	FloatType g,h,fabsh,fabsOffDi,t,theta,c,s,tau,ta,OffDq,a,b;
	const int nxt[] = {Y,Z,X};
	int sweep, i, j;
	U = Matrix4::Identity();
	Diag[X] = S(X,X); Diag[Y] = S(Y,Y); Diag[Z] = S(Z,Z);
	OffD[X] = S(Y,Z); OffD[Y] = S(Z,X); OffD[Z] = S(X,Y);
	for (sweep=20; sweep>0; sweep--) {
		FloatType sm = abs(OffD[X])+abs(OffD[Y])+abs(OffD[Z]);
		if (sm==0.0) break;
		for (i=Z; i>=X; i--) {
			int p = nxt[i]; int q = nxt[p];
			fabsOffDi = abs(OffD[i]);
			g = 100.0*fabsOffDi;
			if (fabsOffDi>0.0) {
				h = Diag[q] - Diag[p];
				fabsh = abs(h);
				if (fabsh+g==fabsh) {
					t = OffD[i]/h;
				} else {
					theta = 0.5*h/OffD[i];
					t = 1.0/(abs(theta)+sqrt(theta*theta+1.0));
					if (theta<0.0) t = -t;
				}
				c = 1.0/sqrt(t*t+1.0); s = t*c;
				tau = s/(c+1.0);
				ta = t*OffD[i]; OffD[i] = 0.0;
				Diag[p] -= ta; Diag[q] += ta;
				OffDq = OffD[q];
				OffD[q] -= s*(OffD[p] + tau*OffD[q]);
				OffD[p] += s*(OffDq   - tau*OffD[p]);
				for (j=Z; j>=X; j--) {
					a = U(j,p); b = U(j,q);
					U(j,p) -= s*(b + tau*a);
					U(j,q) += s*(a - tau*b);
				}
			}
		}
	}
	return Vector3(Diag[X], Diag[Y], Diag[Z]);
}

/******* Spectral Axis Adjustment *******/

/* Given a unit quaternion, q, and a scale vector, k, find a unit quaternion, p,
 * which permutes the axes and turns freely in the plane of duplicate scale
 * factors, such that q p has the largest possible w component, i.e. the
 * smallest possible angle. Permutes k's components to go with q p instead of q.
 * See Ken Shoemake and Tom Duff. AffineTransformation Animation and Polar Decomposition.
 * Proceedings of Graphics Interface 1992. Details on p. 262-263.
 */
Quaternion snuggle(Quaternion q, Vector3& k)
{
#define SQRTHALF ((FloatType)0.7071067811865475244)
#define sgn(n,v)    ((n)?-(v):(v))
#define swap(a,i,j) {a[3]=a[i]; a[i]=a[j]; a[j]=a[3];}
#define cycle(a,p)  if (p) {a[3]=a[0]; a[0]=a[1]; a[1]=a[2]; a[2]=a[3];} else   {a[3]=a[2]; a[2]=a[1]; a[1]=a[0]; a[0]=a[3];}
    Quaternion p;
    FloatType ka[4];
    int i, turn = -1;
    ka[X] = k.x(); ka[Y] = k.y(); ka[Z] = k.z();
    if (ka[X]==ka[Y]) {if (ka[X]==ka[Z]) turn = W; else turn = Z;}
    else {if (ka[X]==ka[Z]) turn = Y; else if (ka[Y]==ka[Z]) turn = X;}
    if (turn>=0) {
	Quaternion qtoz, qp;
	unsigned neg[3], win;
	FloatType mag[3], t;
	static Quaternion qxtoz(0,SQRTHALF,0,SQRTHALF);
	static Quaternion qytoz(SQRTHALF,0,0,SQRTHALF);
	static Quaternion qppmm( 0.5, 0.5,-0.5,-0.5);
	static Quaternion qpppp( 0.5, 0.5, 0.5, 0.5);
	static Quaternion qmpmm(-0.5, 0.5,-0.5,-0.5);
	static Quaternion qpppm( 0.5, 0.5, 0.5,-0.5);
	static Quaternion q0001( 0.0, 0.0, 0.0, 1.0);
	static Quaternion q1000( 1.0, 0.0, 0.0, 0.0);
	switch (turn) {
	default: return (Qt_Conj(q));
	case X: q = q * (qtoz = qxtoz); swap(ka,X,Z) break;
	case Y: q = q * (qtoz = qytoz); swap(ka,Y,Z) break;
	case Z: qtoz = q0001; break;
	}
	q = Qt_Conj(q);
	mag[0] = (FloatType)q.z()*q.z()+(FloatType)q.w()*q.w()-0.5;
	mag[1] = (FloatType)q.x()*q.z()-(FloatType)q.y()*q.w();
	mag[2] = (FloatType)q.y()*q.z()+(FloatType)q.x()*q.w();
	for (i=0; i<3; i++) if ((neg[i] = (mag[i]<0.0))) mag[i] = -mag[i];
	if (mag[0]>mag[1]) {if (mag[0]>mag[2]) win = 0; else win = 2;}
	else		   {if (mag[1]>mag[2]) win = 1; else win = 2;}
	switch (win) {
	case 0: if (neg[0]) p = q1000; else p = q0001; break;
	case 1: if (neg[1]) p = qppmm; else p = qpppp; cycle(ka,0) break;
	case 2: if (neg[2]) p = qmpmm; else p = qpppm; cycle(ka,1) break;
	}
	qp = q * p;
	t = sqrt(mag[win]+0.5);
	p = p * Quaternion(0.0,0.0,-qp.z()/t,qp.w()/t);
	p = qtoz * p.inverse();
    } else {
	FloatType qa[4], pa[4];
	unsigned lo, hi, neg[4], par = 0;
	FloatType all, big, two;
	qa[0] = q.x(); qa[1] = q.y(); qa[2] = q.z(); qa[3] = q.w();
	for (i=0; i<4; i++) {
	    pa[i] = 0.0;
	    if ((neg[i] = (qa[i]<0.0))) qa[i] = -qa[i];
	    par ^= neg[i];
	}
	/* Find two largest components, indices in hi and lo */
	if (qa[0]>qa[1]) lo = 0; else lo = 1;
	if (qa[2]>qa[3]) hi = 2; else hi = 3;
	if (qa[lo]>qa[hi]) {
	    if (qa[lo^1]>qa[hi]) {hi = lo; lo ^= 1;}
	    else {hi ^= lo; lo ^= hi; hi ^= lo;}
	} else {if (qa[hi^1]>qa[lo]) lo = hi^1;}
	all = (qa[0]+qa[1]+qa[2]+qa[3])*0.5;
	two = (qa[hi]+qa[lo])*SQRTHALF;
	big = qa[hi];
	if (all>two) {
	    if (all>big) {/*all*/
		{int i; for (i=0; i<4; i++) pa[i] = sgn(neg[i], 0.5);}
		cycle(ka,par)
	    } else {/*big*/ pa[hi] = sgn(neg[hi],1.0);}
	} else {
	    if (two>big) {/*two*/
		pa[hi] = sgn(neg[hi],SQRTHALF); pa[lo] = sgn(neg[lo], SQRTHALF);
		if (lo>hi) {hi ^= lo; lo ^= hi; hi ^= lo;}
		if (hi==W) {hi = "\001\002\000"[lo]; lo = 3-hi-lo;}
		swap(ka,hi,lo)
	    } else {/*big*/ pa[hi] = sgn(neg[hi],1.0);}
	}
	p.x() = -pa[0]; p.y() = -pa[1]; p.z() = -pa[2]; p.w() = pa[3];
    }
    k.x() = ka[X]; k.y() = ka[Y]; k.z() = ka[Z];
    return (p);
}

/******* Decompose Affine AffineTransformation *******/

/* Decompose 4x4 affine matrix A as TFRUK(U transpose), where t contains the
 * translation components, q contains the rotation R, u contains U, k contains
 * scale factors, and f contains the sign of the determinant.
 * Assumes A transforms column vectors in right-handed coordinates.
 * See Ken Shoemake and Tom Duff. AffineTransformation Animation and Polar Decomposition.
 * Proceedings of Graphics Interface 1992.
 */
void decomp_affine(Matrix4& A, AffineDecomposition* parts)
{
    Matrix4 Q, S, U;
    Quaternion p;
    FloatType det;
    parts->translation = Vector3(A(X,W), A(Y,W), A(Z,W));
    det = polar_decomp(A, Q, S);
    if(det < 0.0) {
		mat_copy(Q,=,-Q,3);
		parts->sign = -1;
    } 
	else parts->sign = 1;
    parts->rotation = Qt_FromMatrix(Q);
	parts->scaling.S = spect_decomp(S, U);
	parts->scaling.Q = Qt_FromMatrix(U);
    p = snuggle(parts->scaling.Q, parts->scaling.S);
    parts->scaling.Q = (parts->scaling.Q * p).normalized();
}

};
