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

import de.jstacs.classifiers.neuralNetwork.activationFunctions.ActivationFunction;
import de.jstacs.data.sequences.Sequence;


public class MSEOutputNeuron extends InnerNeuron {

	/**
	 * @param activationFunction
	 * @param index
	 * @param predecessors
	 */
	public MSEOutputNeuron( ActivationFunction activationFunction, int index, Neuron... predecessors ) {
		super( activationFunction, index, predecessors );
	}
	

	@Override
	public double getError( Sequence input, double weight, double[] desiredOutputs ) {
		if(this.last != input){
			throw new RuntimeException();
		}
		if(error == null){
			error = weight*(desiredOutputs[index]-output) * activationFunction.getDerivation( this.h );
		}
		for(int i=0;i<predecessors.length;i++){
			delta[i] += error*predecessors[i].getOutput( input );
		}
		delta[delta.length-1] += error*(-1);
		return error;
	}

}
