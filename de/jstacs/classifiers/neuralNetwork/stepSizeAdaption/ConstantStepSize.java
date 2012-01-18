package de.jstacs.classifiers.neuralNetwork.stepSizeAdaption;

import de.jstacs.Singleton;


public class ConstantStepSize implements StepSizeAdaption, Singleton {
	
	public static final ConstantStepSize SINGLETON = new ConstantStepSize();
	
	private ConstantStepSize(){}
	
	@Override
	public double getStepSize( double initialStepSize, double lastStepSize, int iteration, int epoch ) {
		return initialStepSize;
	}

	@Override
	public StringBuffer toXML() {
		return new StringBuffer();
	}
	
	

}
