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
package de.jstacs.classifiers.neuralNetwork.neurons;

import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;


public class InputNeuron extends Neuron {

	private int position;
	
	public InputNeuron(int position){
		super();
		this.position = position;
	}
	
	@Override
	public int getNumberOfWeights() {
		return 0;
	}

	@Override
	public double getOutput( Sequence input ) {
		return input.continuousVal( position );
	}

	@Override
	public void initializeRandomly() {		
	}

	@Override
	public void reset() {
	}

	@Override
	public double getError( Sequence input, double weight, double[] desiredOutputs ) {
		double err = 0;
		for(int i=0;i<getNumberOfDescendants();i++){
			err += getDescendant( i ).getError( input, weight, desiredOutputs );
		}
		return err;
	}

	@Override
	public void adaptWeights( double eta ) {
	}



	@Override
	protected void extractFurtherInformation( StringBuffer xml ) throws NonParsableException {
		position = XMLParser.extractObjectForTags( xml, "position", int.class );
	}



	@Override
	protected StringBuffer getFurtherInformation( ) {
		StringBuffer xml = new StringBuffer();
		XMLParser.appendObjectWithTags( xml, position, "position" );
		return xml;
	}

	

}
