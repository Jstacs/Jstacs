package de.jstacs.algorithms.optimization;

/**
 * This interface defines methods for functions that are multi-threaded.
 * 
 * @author Jens Keilwagen
 */
public interface MultiThreadedFunction extends Function {

	/**
	 * This method can and should be used to stop all threads if they are not needed any longer.
	 */
	void stopThreads();

	/**
	 * Returns the number of used threads for evaluating the function and for determining the gradient of the function.
	 * 
	 * @return the number of used threads for evaluating the function and for determining the gradient of the function
	 */
	int getNumberOfThreads();

}