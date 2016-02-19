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
	 * @param curr the current value of progress
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
	
	/**
	 * Increases the internal current value.
	 * @param val the value to be added to the current value
	 */
	public void add( double val ) {
		curr+=val;
	}	
}
