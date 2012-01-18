package de.jstacs.classifiers.neuralNetwork.stepSizeAdaption;

import de.jstacs.Storable;


public interface StepSizeAdaption extends Storable {

	public double getStepSize(double initialStepSize, double lastStepSize, int iteration, int epoch);
	
}
