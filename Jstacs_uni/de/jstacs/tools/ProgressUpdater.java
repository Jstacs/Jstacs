package de.jstacs.tools;

/**
 * Interface for a progress bar (or similar) indicating a tool's progress.
 * 
 * @author Jan Grau
 *
 */
public class ProgressUpdater {

	private double last;
	private double curr; 
	
	/**
	 * Creates a new {@link ProgressUpdater}.
	 */
	public ProgressUpdater(){
		this.last = -1;
	}
	
	/**
	 * Sets the value that is reached upon completion of the monitored task.
	 * @param last the final value
	 */
	public void setLast(double last){
		this.last = last;
	}
	
	/**
	 * Sets the value corresponding to the current progress
	 * @param curr
	 */
	public void setCurrent(double curr){
		this.curr = curr;
	}
	
	/**
	 * Returns the percentage of a tool's work that has been completed so far.
	 * @return the percentage
	 */
	public double getPercentage(){
		if(last > 0){
			return curr/last;
		}else{
			return -1;
		}
	}
	
	/**
	 * Sets the progress to indeterminate.
	 */
	public void setIndeterminate(){
		this.last = -1;
	}
	
	/**
	 * Checks if this progress is indeterminate.
	 * @return if this progress is indeterminate
	 */
	public boolean isIndeterminate(){
		return last<0;
	}
	
}
