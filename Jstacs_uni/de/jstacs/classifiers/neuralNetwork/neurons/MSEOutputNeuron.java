package de.jstacs.classifiers.neuralNetwork.neurons;

import de.jstacs.classifiers.neuralNetwork.activationFunctions.ActivationFunction;
import de.jstacs.data.sequences.Sequence;


public class MSEOutputNeuron extends InnerNeuron {

	/**
	 * @param activationFunction
	 * @param outputLayer
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
