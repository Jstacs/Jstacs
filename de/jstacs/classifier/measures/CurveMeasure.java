package de.jstacs.classifier.measures;

import de.jstacs.DataType;
import de.jstacs.NonParsableException;
import de.jstacs.parameters.SimpleParameter;


public abstract class CurveMeasure extends TwoClassAbstractMeasure {

	public CurveMeasure() {
	}

	public CurveMeasure(boolean getCurve) throws Exception{
		loadParameters();
		parameters.get( 0 ).setValue( getCurve );
	}
	
	public CurveMeasure( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}

	@Override
	protected void loadParameters() throws Exception {
		initParameterList();
		parameters.add( new SimpleParameter( DataType.BOOLEAN, "Curve", "Return the curve of the "+getName(), true ) );
	}

}
