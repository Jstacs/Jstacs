package de.jstacs.classifiers.neuralNetwork.activationFunctions;

import de.jstacs.Storable;


public interface ActivationFunction extends Storable{

	public String getName();
	
	public double getValue(double x);
	
	public double getDerivation(double x);
	
	public double getPositiveValue();
	
	public double getNegativeValue();
	
}
