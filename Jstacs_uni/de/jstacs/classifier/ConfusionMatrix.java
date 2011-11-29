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

package de.jstacs.classifier;

import de.jstacs.NonParsableException;
import de.jstacs.Storable;
import de.jstacs.io.XMLParser;

/**
 * This class holds the confusion matrix of a classifier.
 * 
 * @author Jens Keilwagen
 */
public class ConfusionMatrix implements Storable {

	/**
	 * The confusion matrix.
	 */
	private int[][] matrix;

	private int all;

	private static final String XML_TAG = "confusion matrix";

	/**
	 * Creates a new {@link ConfusionMatrix} with a given number of classes.
	 * 
	 * @param classes
	 *            the number of classes
	 */
	public ConfusionMatrix( int classes ) {
		matrix = new int[classes][classes];
		all = 0;
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link ConfusionMatrix} out of its XML representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link ConfusionMatrix} could not be reconstructed out
	 *             of the XML representation (the {@link StringBuffer}
	 *             <code>representation</code> could not be parsed)
	 * 
	 * @see de.jstacs.Storable
	 */
	public ConfusionMatrix( StringBuffer representation ) throws NonParsableException {
		StringBuffer xml = XMLParser.extractForTag( representation, XML_TAG );
		all = XMLParser.extractObjectForTags( xml, "all", int.class );// TODO XMLP14CONV This and (possibly) the following lines have been converted automatically
		matrix = XMLParser.extractObjectForTags( xml, "matrix", int[][].class );
	}

	/**
	 * This method updates the confusion matrix.
	 * 
	 * @param realClass
	 *            the real class index
	 * @param predictedClass
	 *            the predicted class index
	 */
	public void add( int realClass, int predictedClass ) {
		matrix[predictedClass][realClass]++;
		all++;
	}

	/**
	 * This method returns the confusion matrix as a two dimensional
	 * <code>int</code>-array.
	 * 
	 * @return the confusion matrix
	 */
	public int[][] getMatrix() {
		int[][] res = new int[matrix.length][];
		for( int i = 0; i < matrix.length; i++ ) {
			res[i] = matrix[i].clone();
		}
		return res;
	}

	/**
	 * This method returns the classification rate.
	 * 
	 * @return the classification rate
	 */
	public double getClassificationRate() {
		double correct = 0;
		for( int i = 0; i < matrix.length; i++ ) {
			correct += matrix[i][i];
		}
		return correct / (double) all;
	}
	
	/**
	 * This method returns the misclassification rate.
	 * 
	 * @return the misclassification rate
	 */
	public double getMisclassificationRate() {
		return 1d - getClassificationRate();
	}
	
	/**
	 * This method returns the specific entry of the {@link ConfusionMatrix}.
	 * 
	 * @param predictedClass the index of the predicted class
	 * @param realClass the index of the real class
	 * 
	 * @return the specific entry of the {@link ConfusionMatrix}
	 */
	public int getCountsFor( int predictedClass, int realClass ) {
		return matrix[predictedClass][realClass];
	}

	/* (non-Javadoc)
	 * @see de.jstacs.Storable#toXML()
	 */
	public StringBuffer toXML() {
		StringBuffer xml = new StringBuffer( 1000 );
		XMLParser.appendObjectWithTags( xml, all, "all" );
		XMLParser.appendObjectWithTags( xml, matrix, "matrix" );
		XMLParser.addTags( xml, XML_TAG );
		return xml;
	}
}
