/*
 * This file is part of Jstacs.
 *
 * Jstacs is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Jstacs is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Jstacs.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package de.jstacs.algorithms.optimization.termination;

import de.jstacs.InstantiableFromParameterSet;
import de.jstacs.Storable;
import de.jstacs.utils.Time;

/**
 * This interface can be used in any iterative algorithm for determining the end of the algorithm. 
 * For this reason, the interface declares the method {@link #doNextIteration(int, double, double, double[], double[], double, Time)}
 * which method returns <code>false</code> if no further iteration of the algorithm should be computed and the algorithm should be stopped.
 * If the method returns <code>true</code> another iteration in the algorithm should be done. 
 * 
 * @author Jens Keilwagen
 * 
 * @see de.jstacs.algorithms.optimization.Optimizer
 */
public interface TerminationCondition extends Cloneable, InstantiableFromParameterSet, Storable {
	
	/**
	 * This method allows to decide whether to do another iteration in an optimization or not.
	 * If it returns <code>true</code> it is recommended to do another iteration, otherwise
	 * (the method returns <code>false</code> and) it is recommended to stop the algorithm.
	 * 
	 * @param iteration the number of performed iterations
	 * @param f_last last value of the function
	 * @param f_current current value of the function 
	 * @param gradient the gradient of the function
	 * @param direction the last direction of the optimization
	 * @param alpha the last step size
	 * @param t a time object measuring the time that has been elapsed in the optimization
	 * 
	 * @return <code>true</code> if another iteration should be done
	 */
	boolean doNextIteration( int iteration, double f_last, double f_current, double[] gradient, double[] direction, double alpha, Time t );
	
	/**
	 * This method returns <code>false</code> if the {@link TerminationCondition} uses either
	 * the gradient or the direction for the decision, otherwise it returns <code>true</code>.
	 * 
	 * @return <code>false</code> for gradient or direction based {@link TerminationCondition}s
	 */
	boolean isSimple();
}
