package de.jstacs.clustering.hierachical;


/**
 * Interface for classes that may provide a position weight matrix
 * 
 * @author Jan Grau
 *
 */
public interface PWMSupplier {

	/**
	 * Returns the position weight matrix. Rows are positions, columns are symbols 
	 * @return the PWM
	 */
	public double[][] getPWM();
	
	/**
	 * Returns a name (e.g., an identifier from a database) for the PWM.
	 * @return the name
	 */
	public String getName();
	
}
