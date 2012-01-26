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
