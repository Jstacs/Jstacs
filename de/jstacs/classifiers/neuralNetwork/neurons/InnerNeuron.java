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

import java.util.Arrays;
import java.util.Random;

import de.jstacs.classifiers.neuralNetwork.activationFunctions.ActivationFunction;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;


public class InnerNeuron extends Neuron {

	private static Random r = new Random();
	
	protected Neuron[] predecessors;
	protected double[] weights;
	protected int index;
	protected ActivationFunction activationFunction;
	protected Double h = null;
	protected Double output = null;
	protected Double error = null;
	protected double[] delta;
	protected Sequence last;
	
	public InnerNeuron(ActivationFunction activationFunction, int index, Neuron... predecessors){
		super();
		this.predecessors = predecessors.clone();//do not clone instances
		for(int i=0;i<this.predecessors.length;i++){
			predecessors[i].addDescendant( this );
		}
		this.weights = new double[predecessors.length+1];
		this.delta = new double[weights.length];
		this.activationFunction = activationFunction;
		this.index = index;
	}
	
	public final void setPredecessors(Neuron... predecessors){//TODO security
		if(this.predecessors != null){
			throw new IllegalArgumentException( "Can set predecessors only initially" );
		}
		this.predecessors = predecessors.clone();
		for(int i=0;i<this.predecessors.length;i++){
			predecessors[i].addDescendant( this );
		}
	}
	
	@Override
	public int getNumberOfWeights() {
		return weights.length;
	}
	
	@Override
	public double getOutput( Sequence input ) {
		if(last != null && input != last){
			reset();
		}
		if(this.output == null){
			double output = 0.0;
			for(int i=0;i<predecessors.length;i++){
				output += weights[i]*predecessors[i].getOutput( input );
			}
			output -= weights[weights.length-1];
			this.h = output;
			this.output = activationFunction.getValue( output );
			this.last = input;
		}
		return this.output;
	}

	public void reset(){
		output = null;
		error = null;
		h = null;
	}

	@Override
	public void initializeRandomly() {
		for(int i=0;i<weights.length;i++){
			weights[i] = r.nextGaussian();
		}
	}

	@Override
	public double getError( Sequence input, double weight, double[] desiredOutputs ) {
		if(this.last != input){
			throw new RuntimeException();
		}
		if(error == null){
			double err = 0;
			for(int i=0;i<getNumberOfDescendants();i++){
				err += getDescendant( i ).getError( input, weight, desiredOutputs )*getDescendant( i ).getWeightForPredecessor( index );
			}
			error = err * activationFunction.getDerivation( this.h );
			for(int i=0;i<predecessors.length;i++){
				delta[i] += error*predecessors[i].getOutput( input );
			}
			delta[delta.length-1] += error*(-1);
		}
		return error;
	}

	@Override
	public void adaptWeights(double eta){
		
		for(int i=0;i<weights.length;i++){
			weights[i] += eta*delta[i];
		}
		Arrays.fill( delta, 0 );
	}

	public double getWeightForPredecessor( int index ) {
		return weights[index];
	}

	@Override
	protected void extractFurtherInformation( StringBuffer xml ) throws NonParsableException {
		activationFunction = XMLParser.extractObjectForTags( xml, "activationFunction", ActivationFunction.class );
		delta = XMLParser.extractObjectForTags( xml, "delta", double[].class );
		index = XMLParser.extractObjectForTags( xml, "index", int.class );
		weights = XMLParser.extractObjectForTags( xml, "weights", double[].class );
	}

	@Override
	protected StringBuffer getFurtherInformation( ) {
		StringBuffer xml = new StringBuffer();
		XMLParser.appendObjectWithTags( xml, activationFunction, "activationFunction" );
		XMLParser.appendObjectWithTags( xml, delta, "delta" );
		XMLParser.appendObjectWithTags( xml, index, "index" );
		XMLParser.appendObjectWithTags( xml, weights, "weights" );
		return xml;
	}	
	
	
	
}
