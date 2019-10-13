////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright 2018 Alexander Stukowski
//
//  This file is part of OVITO (Open Visualization Tool).
//
//  OVITO is free software; you can redistribute it and/or modify it either under the
//  terms of the GNU General Public License version 3 as published by the Free Software
//  Foundation (the "GPL") or, at your option, under the terms of the MIT License.
//  If you do not alter this notice, a recipient may use your version of this
//  file under either the GPL or the MIT License.
//
//  You should have received a copy of the GPL along with this program in a
//  file LICENSE.GPL.txt.  You should have received a copy of the MIT License along
//  with this program in a file LICENSE.MIT.txt
//
//  This software is distributed on an "AS IS" basis, WITHOUT WARRANTY OF ANY KIND,
//  either express or implied. See the GPL or the MIT License for the specific language
//  governing rights and limitations.
//
////////////////////////////////////////////////////////////////////////////////////////

#pragma once


#include <ovito/core/Core.h>
#include <ovito/core/oo/RefTarget.h>
#include <ovito/core/dataset/pipeline/PipelineStatus.h>
#include <ovito/core/dataset/animation/TimeInterval.h>

namespace Ovito { OVITO_BEGIN_INLINE_NAMESPACE(View)

/**
 * \brief Abstract base class for all viewport overlays.
 */
class OVITO_CORE_EXPORT ViewportOverlay : public RefTarget
{
	Q_OBJECT
	OVITO_CLASS(ViewportOverlay)

protected:

	/// \brief Constructor.
	ViewportOverlay(DataSet* dataset);

public:

	/// \brief This method asks the overlay to paint its contents over the rendered image.
	virtual void render(const Viewport* viewport, TimePoint time, FrameBuffer* frameBuffer,
						const ViewProjectionParameters& projParams, const RenderSettings* renderSettings, AsyncOperation& operation) = 0;

	/// \brief This method asks the overlay to paint its contents over the given interactive viewport.
	virtual void renderInteractive(const Viewport* viewport, TimePoint time, QPainter& painter,
						const ViewProjectionParameters& projParams, const RenderSettings* renderSettings, AsyncOperation& operation) = 0;


	/// \brief Lets the overlay render its 3d content.
	/// \param vp The viewport into which to render the graphical content.
	/// \param renderer The renderer that should be used to produce the visualization.
	///
	/// The default implementation of this method does nothing.
	virtual void render3D(Viewport* vp, TimePoint time, SceneRenderer* renderer, AsyncOperation& operation) {}

	/// \brief Moves the position of the overlay in the viewport by the given amount,
	///        which is specified as a fraction of the viewport render size.
	///
	/// Overlay implementations should override this method if they support positioning.
	/// The default method implementation does nothing.
	virtual void moveOverlayInViewport(const Vector2& delta) {}

protected:

	/// \brief Is called when the value of a non-animatable property field of this RefMaker has changed.
	virtual void propertyChanged(const PropertyFieldDescriptor& field) override;

private:

	/// The current status of this overlay object.
	DECLARE_RUNTIME_PROPERTY_FIELD_FLAGS(PipelineStatus, status, setStatus, PROPERTY_FIELD_NO_UNDO | PROPERTY_FIELD_NO_CHANGE_MESSAGE);

	/// Option for rendering the overlay contents behind the three-dimensional content.
	DECLARE_MODIFIABLE_PROPERTY_FIELD(bool, renderBehindScene, setRenderBehindScene);

	/// Flag controlling the visibility of the overlay.
	DECLARE_MODIFIABLE_PROPERTY_FIELD(bool, isEnabled, setEnabled);
};

OVITO_END_INLINE_NAMESPACE
}	// End of namespace
