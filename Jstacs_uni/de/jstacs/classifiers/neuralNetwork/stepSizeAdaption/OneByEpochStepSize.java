package de.jstacs.classifiers.neuralNetwork.stepSizeAdaption;

import de.jstacs.Singleton;


public class OneByEpochStepSize implements StepSizeAdaption, Singleton {

	public static final OneByEpochStepSize SINGLETON = new OneByEpochStepSize();
	
	private OneByEpochStepSize(){}
	
	@Override
	public double getStepSize( double initialStepSize, double lastStepSize, int iteration, int epoch ) {
		return initialStepSize/(double)epoch;
	}

	@Override
	public StringBuffer toXML() {
		return new StringBuffer();
	}
	
}
