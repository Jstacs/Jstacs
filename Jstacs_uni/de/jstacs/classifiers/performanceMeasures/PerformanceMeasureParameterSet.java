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
 * This class implements a container of {@link AbstractPerformanceMeasure}s that can be used
 * in {@link de.jstacs.classifiers.AbstractClassifier#evaluate(PerformanceMeasureParameterSet, boolean, de.jstacs.data.DataSet...)}.
 * 
 * @author Jens Keilwagen, Jan Grau
 */
public class PerformanceMeasureParameterSet extends AbstractPerformanceMeasureParameterSet<PerformanceMeasure> {

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link PerformanceMeasureParameterSet} out of an XML representation.
	 * 
	 * @param xml
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link PerformanceMeasureParameterSet} could not be reconstructed out of
	 *             the {@link StringBuffer} <code>xml</code>
	 */
	public PerformanceMeasureParameterSet( StringBuffer xml ) throws NonParsableException {
		super( xml );
	}
		
	/**
	 * Constructs a new {@link PerformanceMeasureParameterSet} that can be used for classifiers that
	 * handle the given number of classes.
	 * 
	 * @param numClasses the number of classes
	 *  
	 * @throws Exception if something went wrong
	 * 
	 * @see de.jstacs.classifiers.AbstractClassifier#getNumberOfClasses()
	 */
	public PerformanceMeasureParameterSet( int numClasses ) throws Exception {
		super( numClasses, false, new PerformanceMeasure[0] );
	}
	
	/**
	 * Constructs a new {@link PerformanceMeasureParameterSet} with the given performance measures.
	 * The number of classes this {@link PerformanceMeasureParameterSet} can be used for is determined from
	 * the given {@link AbstractPerformanceMeasure}s.  If no measure is given, a new
	 * {@link PerformanceMeasureParameterSet} for binary classifiers is created.
	 * 
	 * @param measures the {@link AbstractPerformanceMeasure} that shall be used
	 *  
	 * @throws Exception if something went wrong
	 */
	public PerformanceMeasureParameterSet( PerformanceMeasure... measures ) throws Exception {
		super( getNumberOfClasses( measures ), false, measures );
	}	
}
