/*
 * This file is part of Jstacs.
 *
 * Jstacs is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * Jstacs is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Jstacs. If not, see <http://www.gnu.org/licenses/>.
 *
 * For more information on Jstacs, visit http://www.jstacs.de
 */
package de.jstacs.classifiers.neuralNetwork.activationFunctions;

import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;


public class LogisticFunctionActivationFunction implements ActivationFunction {

	private double beta;
	
	public LogisticFunctionActivationFunction(double beta){
		this.beta = beta;
	}
	
	public LogisticFunctionActivationFunction(StringBuffer xml) throws NonParsableException{
		xml = XMLParser.extractForTag( xml, getClass().getSimpleName() );
		beta = XMLParser.extractObjectForTags( xml, "beta", double.class );
	}
	
	@Override
	public String getName() {
		return "Logistic function with beta="+beta;
	}

	@Override
	public double getValue( double x ) {
		return 1.0/(1.0+Math.exp( -beta*x ));
	}

	@Override
	public double getDerivation( double x ) {
		double a = getValue( x );
		return beta*a*(1.0-a);
	}

	@Override
	public double getPositiveValue() {
		return 1;
	}

	@Override
	public double getNegativeValue() {
		return 0;
	}

	public StringBuffer toXML(){
		StringBuffer xml = new StringBuffer();
		XMLParser.appendObjectWithTags( xml, beta, "beta" );
		XMLParser.addTags( xml, getClass().getSimpleName() );
		return xml;
	}
	
}
