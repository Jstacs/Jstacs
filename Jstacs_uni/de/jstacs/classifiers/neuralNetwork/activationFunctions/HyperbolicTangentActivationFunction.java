package de.jstacs.classifiers.neuralNetwork.activationFunctions;

import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;


public class HyperbolicTangentActivationFunction implements ActivationFunction {

	private double beta;
	
	public HyperbolicTangentActivationFunction(double beta){
		this.beta = beta;
	}
	
	public HyperbolicTangentActivationFunction(StringBuffer xml) throws NonParsableException{
		xml = XMLParser.extractForTag( xml, getClass().getSimpleName() );
		beta = XMLParser.extractObjectForTags( xml, "beta", double.class );
	}
	
	@Override
	public String getName() {
		return "Hyperbolic tangent with beta="+beta;
	}

	@Override
	public double getValue( double x ) {
		return Math.tanh( beta*x );
	}

	@Override
	public double getDerivation( double x ) {
		double a = Math.tanh( beta*x );
		return beta*(1.0-a*a);
	}

	@Override
	public double getPositiveValue() {
		return 1;
	}

	@Override
	public double getNegativeValue() {
		return -1;
	}
	
	public StringBuffer toXML(){
		StringBuffer xml = new StringBuffer();
		XMLParser.appendObjectWithTags( xml, beta, "beta" );
		XMLParser.addTags( xml, getClass().getSimpleName() );
		return xml;
	}

}
