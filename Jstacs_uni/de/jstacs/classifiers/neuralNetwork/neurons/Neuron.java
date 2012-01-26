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

import java.util.ArrayList;

import de.jstacs.Storable;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;


public abstract class Neuron implements Storable{
	
	private ArrayList<InnerNeuron> descendants;
	
	public Neuron(){
		this.descendants = new ArrayList<InnerNeuron>();
	}
	
	public Neuron(StringBuffer xml) throws NonParsableException{
		fromXML(xml);
		if(this.descendants == null){
			throw new NonParsableException();
		}
	}
	
	private void fromXML( StringBuffer xml ) throws NonParsableException {
		xml = XMLParser.extractForTag( xml, getClass().getSimpleName() );
		extractFurtherInformation(xml);
	}

	protected abstract void extractFurtherInformation( StringBuffer xml ) throws NonParsableException;

	protected abstract StringBuffer getFurtherInformation( );

	public StringBuffer toXML(){
		StringBuffer sb = new StringBuffer();
		sb.append( getFurtherInformation( ));
		XMLParser.addTags( sb, getClass().getSimpleName() );
		return sb;
	}
	
	public final int getNumberOfDescendants(){
		return descendants.size();
	}
	
	public final InnerNeuron getDescendant(int i){
		return descendants.get( i );
	}
	
	protected final void addDescendant( InnerNeuron neuron ){
		this.descendants.add( neuron );
	}
		
	public abstract int getNumberOfWeights();
	
	public abstract double getOutput(Sequence input);
	
	public abstract void initializeRandomly();
	
	public abstract void reset();
	
	public abstract double getError(Sequence input, double weight, double[] desiredOutputs );

	public abstract void adaptWeights( double eta );
}
