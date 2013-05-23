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
package de.jstacs.classifiers.performanceMeasures;

import de.jstacs.io.NonParsableException;

/**
 * This class implements a container for {@link NumericalPerformanceMeasure}s that can be used, for instance, in an repeated assessment,
 * (cf. {@link de.jstacs.classifiers.assessment.ClassifierAssessment}).
 * 
 * @author Jens Keilwagen
 */
public class NumericalPerformanceMeasureParameterSet extends AbstractPerformanceMeasureParameterSet<NumericalPerformanceMeasure> {

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link NumericalPerformanceMeasureParameterSet} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link NumericalPerformanceMeasureParameterSet} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
	public NumericalPerformanceMeasureParameterSet( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}
	
	/**
	 * Constructs a new {@link NumericalPerformanceMeasureParameterSet} that can be used for classifiers that
	 * handle the given number of classes. Automatically includes all {@link AbstractPerformanceMeasure}s that can be
	 * computed for the given number of classes.
	 * 
	 * @param numClasses the number of classes
	 *  
	 * @throws Exception if something went wrong
	 * 
	 * @see de.jstacs.classifiers.AbstractClassifier#getNumberOfClasses()
	 */
	public NumericalPerformanceMeasureParameterSet( int numClasses ) throws Exception {
		super( numClasses, true, new NumericalPerformanceMeasure[0] );
	}
	
	/**
	 * Constructs a new {@link NumericalPerformanceMeasureParameterSet} with the given performance measures.
	 * The number of classes this {@link NumericalPerformanceMeasureParameterSet} can be used for is determined from
	 * the given {@link NumericalPerformanceMeasure}s. If no measure is given, a new
	 * {@link NumericalPerformanceMeasureParameterSet} for binary classifiers is created.
	 * 
	 * @param measures the {@link NumericalPerformanceMeasure} that shall be used
	 *  
	 * @throws Exception if something went wrong
	 */
	public NumericalPerformanceMeasureParameterSet( NumericalPerformanceMeasure... measures ) throws Exception {
		super( getNumberOfClasses( measures ), false, measures );
	}
}
