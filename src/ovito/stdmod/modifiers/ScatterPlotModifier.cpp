///////////////////////////////////////////////////////////////////////////////
//
//  Copyright (2017) Alexander Stukowski
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
///////////////////////////////////////////////////////////////////////////////

#include <ovito/stdmod/StdMod.h>
#include <ovito/core/dataset/DataSet.h>
#include <ovito/core/dataset/pipeline/ModifierApplication.h>
#include <ovito/stdobj/properties/PropertyObject.h>
#include <ovito/stdobj/properties/PropertyContainer.h>
#include <ovito/stdobj/series/DataSeriesObject.h>
#include <ovito/core/app/Application.h>
#include "ScatterPlotModifier.h"

namespace Ovito { namespace StdMod {

IMPLEMENT_OVITO_CLASS(ScatterPlotModifier);
DEFINE_PROPERTY_FIELD(ScatterPlotModifier, selectXAxisInRange);
DEFINE_PROPERTY_FIELD(ScatterPlotModifier, selectionXAxisRangeStart);
DEFINE_PROPERTY_FIELD(ScatterPlotModifier, selectionXAxisRangeEnd);
DEFINE_PROPERTY_FIELD(ScatterPlotModifier, selectYAxisInRange);
DEFINE_PROPERTY_FIELD(ScatterPlotModifier, selectionYAxisRangeStart);
DEFINE_PROPERTY_FIELD(ScatterPlotModifier, selectionYAxisRangeEnd);
DEFINE_PROPERTY_FIELD(ScatterPlotModifier, fixXAxisRange);
DEFINE_PROPERTY_FIELD(ScatterPlotModifier, xAxisRangeStart);
DEFINE_PROPERTY_FIELD(ScatterPlotModifier, xAxisRangeEnd);
DEFINE_PROPERTY_FIELD(ScatterPlotModifier, fixYAxisRange);
DEFINE_PROPERTY_FIELD(ScatterPlotModifier, yAxisRangeStart);
DEFINE_PROPERTY_FIELD(ScatterPlotModifier, yAxisRangeEnd);
DEFINE_PROPERTY_FIELD(ScatterPlotModifier, xAxisProperty);
DEFINE_PROPERTY_FIELD(ScatterPlotModifier, yAxisProperty);
SET_PROPERTY_FIELD_LABEL(ScatterPlotModifier, selectXAxisInRange, "Select elements in x-range");
SET_PROPERTY_FIELD_LABEL(ScatterPlotModifier, selectionXAxisRangeStart, "Selection x-range start");
SET_PROPERTY_FIELD_LABEL(ScatterPlotModifier, selectionXAxisRangeEnd, "Selection x-range end");
SET_PROPERTY_FIELD_LABEL(ScatterPlotModifier, selectYAxisInRange, "Select elements in y-range");
SET_PROPERTY_FIELD_LABEL(ScatterPlotModifier, selectionYAxisRangeStart, "Selection y-range start");
SET_PROPERTY_FIELD_LABEL(ScatterPlotModifier, selectionYAxisRangeEnd, "Selection y-range end");
SET_PROPERTY_FIELD_LABEL(ScatterPlotModifier, fixXAxisRange, "Fix x-range");
SET_PROPERTY_FIELD_LABEL(ScatterPlotModifier, xAxisRangeStart, "X-range start");
SET_PROPERTY_FIELD_LABEL(ScatterPlotModifier, xAxisRangeEnd, "X-range end");
SET_PROPERTY_FIELD_LABEL(ScatterPlotModifier, fixYAxisRange, "Fix y-range");
SET_PROPERTY_FIELD_LABEL(ScatterPlotModifier, yAxisRangeStart, "Y-range start");
SET_PROPERTY_FIELD_LABEL(ScatterPlotModifier, yAxisRangeEnd, "Y-range end");
SET_PROPERTY_FIELD_LABEL(ScatterPlotModifier, xAxisProperty, "X-axis property");
SET_PROPERTY_FIELD_LABEL(ScatterPlotModifier, yAxisProperty, "Y-axis property");

/******************************************************************************
* Constructs the modifier object.
******************************************************************************/
ScatterPlotModifier::ScatterPlotModifier(DataSet* dataset) : GenericPropertyModifier(dataset),
	_selectXAxisInRange(false),
	_selectionXAxisRangeStart(0),
	_selectionXAxisRangeEnd(1),
	_selectYAxisInRange(false),
	_selectionYAxisRangeStart(0),
	_selectionYAxisRangeEnd(1),
	_fixXAxisRange(false),
	_xAxisRangeStart(0),
	_xAxisRangeEnd(0),
	_fixYAxisRange(false),
	_yAxisRangeStart(0),
	_yAxisRangeEnd(0)
{
	// Operate on particle properties by default.
	setDefaultSubject(QStringLiteral("Particles"), QStringLiteral("ParticlesObject"));
}

/******************************************************************************
* This method is called by the system when the modifier has been inserted
* into a pipeline.
******************************************************************************/
void ScatterPlotModifier::initializeModifier(ModifierApplication* modApp)
{
	GenericPropertyModifier::initializeModifier(modApp);

	// Use the first available property from the input state as data source when the modifier is newly created.
	if((xAxisProperty().isNull() || yAxisProperty().isNull()) && subject() && Application::instance()->executionContext() == Application::ExecutionContext::Interactive) {
		const PipelineFlowState& input = modApp->evaluateInputPreliminary();
		if(const PropertyContainer* container = input.getLeafObject(subject())) {
			PropertyReference bestProperty;
			for(PropertyObject* property : container->properties()) {
				bestProperty = PropertyReference(subject().dataClass(), property, (property->componentCount() > 1) ? 0 : -1);
			}
			if(xAxisProperty().isNull() && !bestProperty.isNull()) {
				setXAxisProperty(bestProperty);
			}
			if(yAxisProperty().isNull() && !bestProperty.isNull()) {
				setYAxisProperty(bestProperty);
			}
		}
	}
}

/******************************************************************************
* Is called when the value of a property of this object has changed.
******************************************************************************/
void ScatterPlotModifier::propertyChanged(const PropertyFieldDescriptor& field)
{
	// Whenever the selected property class of this modifier is changed, update the source property references.
	if(field == PROPERTY_FIELD(GenericPropertyModifier::subject) && !isBeingLoaded() && !dataset()->undoStack().isUndoingOrRedoing()) {
		setXAxisProperty(xAxisProperty().convertToContainerClass(subject().dataClass()));
		setYAxisProperty(yAxisProperty().convertToContainerClass(subject().dataClass()));
	}
	GenericPropertyModifier::propertyChanged(field);
}

/******************************************************************************
* Modifies the input data in an immediate, preliminary way.
******************************************************************************/
void ScatterPlotModifier::evaluatePreliminary(TimePoint time, ModifierApplication* modApp, PipelineFlowState& state)
{
	if(!subject())
		throwException(tr("No data element type set."));
	if(xAxisProperty().isNull())
		throwException(tr("No input property for x-axis selected."));
	if(yAxisProperty().isNull())
		throwException(tr("No input property for y-axis selected."));

	// Check if the source property is the right kind of property.
	if(xAxisProperty().containerClass() != subject().dataClass())
		throwException(tr("Modifier was set to operate on '%1', but the selected input is a '%2' property.")
			.arg(subject().dataClass()->pythonName()).arg(xAxisProperty().containerClass()->propertyClassDisplayName()));

	// Check if the source property is the right kind of property.
	if(yAxisProperty().containerClass() != subject().dataClass())
		throwException(tr("Modifier was set to operate on '%1', but the selected input is a '%2' property.")
			.arg(subject().dataClass()->pythonName()).arg(yAxisProperty().containerClass()->propertyClassDisplayName()));

	// Look up the property container object.
	const PropertyContainer* container = state.expectLeafObject(subject());

	// Get the input properties.
	const PropertyObject* xPropertyObj = xAxisProperty().findInContainer(container);
	if(!xPropertyObj)
		throwException(tr("The selected input property '%1' is not present.").arg(xAxisProperty().name()));
	const PropertyObject* yPropertyObj = yAxisProperty().findInContainer(container);
	if(!yPropertyObj)
		throwException(tr("The selected input property '%1' is not present.").arg(yAxisProperty().name()));

	ConstPropertyPtr xProperty = xPropertyObj->storage();
	ConstPropertyPtr yProperty = yPropertyObj->storage();
	OVITO_ASSERT(xProperty->size() == yProperty->size());

	size_t xVecComponent = std::max(0, xAxisProperty().vectorComponent());
	size_t xVecComponentCount = xProperty->componentCount();
	size_t yVecComponent = std::max(0, yAxisProperty().vectorComponent());
	size_t yVecComponentCount = yProperty->componentCount();
	if(xVecComponent >= xProperty->componentCount())
		throwException(tr("The selected vector component is out of range. The property '%1' has only %2 components per element.").arg(xProperty->name()).arg(xProperty->componentCount()));
	if(yVecComponent >= yProperty->componentCount())
		throwException(tr("The selected vector component is out of range. The property '%1' has only %2 components per element.").arg(yProperty->name()).arg(yProperty->componentCount()));

	// Get selection ranges.
	FloatType selectionXAxisRangeStart = this->selectionXAxisRangeStart();
	FloatType selectionXAxisRangeEnd = this->selectionXAxisRangeEnd();
	FloatType selectionYAxisRangeStart = this->selectionYAxisRangeStart();
	FloatType selectionYAxisRangeEnd = this->selectionYAxisRangeEnd();
	if(selectionXAxisRangeStart > selectionXAxisRangeEnd)
		std::swap(selectionXAxisRangeStart, selectionXAxisRangeEnd);
	if(selectionYAxisRangeStart > selectionYAxisRangeEnd)
		std::swap(selectionYAxisRangeStart, selectionYAxisRangeEnd);

	// Create output selection property.
	PropertyPtr outputSelection;
	size_t numSelected = 0;
	if(selectXAxisInRange() || selectYAxisInRange()) {
		// First make sure we can safely modify the property container.
		PropertyContainer* mutableContainer = state.expectMutableLeafObject(subject());
		// Add the selection property to the output container.
		outputSelection = mutableContainer->createProperty(PropertyStorage::GenericSelectionProperty)->modifiableStorage();
		std::fill(outputSelection->dataInt(), outputSelection->dataInt() + outputSelection->size(), 1);
		numSelected = outputSelection->size();
	}

#if 0
	FloatType xIntervalStart = xAxisRangeStart();
	FloatType xIntervalEnd = xAxisRangeEnd();
	FloatType yIntervalStart = yAxisRangeStart();
	FloatType yIntervalEnd = yAxisRangeEnd();
#endif

	// Create output arrays.
	auto out_x = DataSeriesObject::OOClass().createStandardStorage(container->elementCount(), DataSeriesObject::XProperty, false);
	auto out_y = DataSeriesObject::OOClass().createStandardStorage(container->elementCount(), DataSeriesObject::YProperty, false);

	// Collect X coordinates.
	if(!xProperty->copyTo(out_x->dataFloat(), xVecComponent))
		throwException(tr("Failed to extract coordinate values from input property for x-axis."));

	// Collect Y coordinates.
	if(!yProperty->copyTo(out_y->dataFloat(), yVecComponent))
		throwException(tr("Failed to extract coordinate values from input property for y-axis."));

#if 0
	// Determine value ranges.
	if(fixXAxisRange() == false || fixYAxisRange() == false) {
		Box2 bbox;
		for(const QPointF& p : xyData) {
			bbox.addPoint(p.x(), p.y());
		}
		if(fixXAxisRange() == false) {
			xIntervalStart = bbox.minc.x();
			xIntervalEnd = bbox.maxc.x();
		}
		if(fixYAxisRange() == false) {
			yIntervalStart = bbox.minc.y();
			yIntervalEnd = bbox.maxc.y();
		}
	}
#endif

	if(outputSelection && selectXAxisInRange()) {
		OVITO_ASSERT(outputSelection->size() == out_x->size());
		int* s = outputSelection->dataInt();
		int* s_end = s + outputSelection->size();
		for(FloatType x : out_x->constFloatRange()) {
			if(x < selectionXAxisRangeStart || x > selectionXAxisRangeEnd) {
				*s = 0;
				numSelected--;
			}
			++s;
		}
	}

	if(outputSelection && selectYAxisInRange()) {
		OVITO_ASSERT(outputSelection->size() == out_y->size());
		int* s = outputSelection->dataInt();
		int* s_end = s + outputSelection->size();
		for(FloatType y : out_y->constFloatRange()) {
			if(y < selectionYAxisRangeStart || y > selectionYAxisRangeEnd) {
				if(*s) {
					*s = 0;
					numSelected--;
				}
			}
			++s;
		}
	}

	// Output a data series object with the scatter points.
	DataSeriesObject* seriesObj = state.createObject<DataSeriesObject>(QStringLiteral("scatter"), modApp, DataSeriesObject::Scatter, tr("%1 vs. %2").arg(yAxisProperty().nameWithComponent()).arg(xAxisProperty().nameWithComponent()));
	seriesObj->createProperty(std::move(out_x))->setName(xAxisProperty().nameWithComponent());
	seriesObj->createProperty(std::move(out_y))->setName(yAxisProperty().nameWithComponent());

	QString statusMessage;
	if(outputSelection) {
		statusMessage = tr("%1 %2 selected (%3%)").arg(numSelected)
				.arg(container->getOOMetaClass().elementDescriptionName())
				.arg((FloatType)numSelected * 100 / std::max((size_t)1,outputSelection->size()), 0, 'f', 1);
	}

	state.setStatus(PipelineStatus(PipelineStatus::Success, std::move(statusMessage)));
}

}	// End of namespace
}	// End of namespace