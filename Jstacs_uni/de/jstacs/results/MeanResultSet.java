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
 * along with Jstacs. If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package de.jstacs.results;

import de.jstacs.AnnotatedEntityList;
import de.jstacs.DataType;
import de.jstacs.Storable;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;
import de.jstacs.utils.DoubleList;

/**
 * Class that computes the mean and the standard error of a series of
 * {@link NumericalResultSet}s. Each {@link NumericalResultSet} in the series
 * must be added using the method {@link #addResults(NumericalResultSet...)}.
 * The means and the standard errors can be obtained by {@link #getStatistics()}
 * .
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public class MeanResultSet extends NumericalResultSet {

	/**
	 * The number of results
	 */
	private int count;

	/**
	 * The squares of each result
	 */
	private double[] squares;

	/**
	 * The information to this {@link MeanResultSet}
	 */
	private SimpleResult[] infos;

	private boolean aggregateAll;
	private DoubleList[] allValues;

	private MeanResultSet(NumericalResult[] res, SimpleResult[] infos,
			double[] squares, int count, boolean aggregateAll, DoubleList[] allValues) {
		super(res);
		this.count = count;
		this.infos = new SimpleResult[infos.length];
		this.squares = squares;
		this.aggregateAll = aggregateAll;
		this.allValues = allValues;		
		System.arraycopy(infos, 0, this.infos, 0, infos.length);
	}

	/**
	 * Constructs a new {@link MeanResultSet} with an empty set of
	 * {@link NumericalResultSet}s and allows to collect all
	 * values via the switch <code>aggregateAll</code>.
	 * 
	 * @param infos
	 *            some information to this {@link MeanResultSet}
	 * @param aggregateAll
	 *			  a switch to decide whether all values should be collected
	 */
	public MeanResultSet(boolean aggregateAll, SimpleResult... infos) {
		super();
		count = 0;
		this.infos = new SimpleResult[infos.length + 1];
		System.arraycopy(infos, 0, this.infos, 0, infos.length);
		this.aggregateAll = aggregateAll;
	}
	
	/**
	 * Constructs a new {@link MeanResultSet} with an empty set of
	 * {@link NumericalResultSet}s.
	 * 
	 * @param infos
	 *            some information to this {@link MeanResultSet}
	 */
	public MeanResultSet(SimpleResult... infos) {
		this( false, infos );
	}

	/**
	 * Constructs a new {@link MeanResultSet} with an empty set of
	 * {@link NumericalResultSet}s and no further information.
	 */
	public MeanResultSet() {
		this(new SimpleResult[0]);
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link MeanResultSet} from the corresponding XML
	 * representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link StringBuffer} <code>representation</code> could
	 *             not be parsed
	 */
	public MeanResultSet(StringBuffer representation)
			throws NonParsableException {
		super(representation);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.results.ResultSet#toXML()
	 */
	@Override
	public StringBuffer toXML() {
		StringBuffer buf = super.toXML();
		XMLParser.addTags(buf, "superSet");
		setCount();
		XMLParser.appendObjectWithTags(buf, infos, "infos");
		XMLParser.appendObjectWithTags( buf, aggregateAll, "aggregateAll" );
		if( count > 0 ) {
			XMLParser.appendObjectWithTags(buf, squares, "squares");
			XMLParser.appendObjectWithTags(buf, allValues, "allValues");
		}
		XMLParser.addTags(buf, "meanResultSet");

		return buf;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.results.ResultSet#fromXML(java.lang.StringBuffer)
	 */
	@Override
	protected void fromXML(StringBuffer representation) throws NonParsableException {
		representation = XMLParser.extractForTag(representation,
				"meanResultSet");
		super.fromXML(XMLParser.extractForTag(representation, "superSet"));
		
		Storable[] infosTemp = XMLParser.extractObjectForTags( representation, "infos", Storable[].class);
		infos = new SimpleResult[infosTemp.length];
		for (int i = 0; i < infos.length; i++) {
			infos[i] = (SimpleResult) infosTemp[i];
		}
		aggregateAll = (Boolean) XMLParser.extractObjectForTags( representation, "aggregateAll" );
		count = (Integer) infos[infos.length - 1].getValue();
		if( count > 0 ) {
			squares = XMLParser.extractObjectForTags(representation, "squares", double[].class );
			allValues = (DoubleList[]) XMLParser.extractObjectForTags( representation, "allValues" );
		}
	}

	/**
	 * Adds two {@link MeanResultSet}s together.
	 * 
	 * @param r1
	 *            one {@link MeanResultSet}
	 * @param r2
	 *            another {@link MeanResultSet}
	 * 
	 * @return the merged {@link MeanResultSet}
	 * 
	 * @throws AdditionImpossibleException
	 *             if the {@link Result}s does not match
	 */
	public static MeanResultSet addResults(MeanResultSet r1, MeanResultSet r2)
			throws AdditionImpossibleException {
		if (r1.getNumberOfResults() != r2.getNumberOfResults()) {
			throw new AdditionImpossibleException();
		}
		int i = 0;
		while (i < r1.infos.length - 1 && r1.infos[i].equals(r2.infos[i])) {
			i++;
		}
		if (i != r1.infos.length - 1) {
			throw new AdditionImpossibleException();
		}
		boolean aggregateAll = r1.aggregateAll && r2.aggregateAll;
		NumericalResult[] results = new NumericalResult[r1.getNumberOfResults()];
		double[] squares = new double[r1.getNumberOfResults()];
		DoubleList[] allValues = null;
		if( aggregateAll ) {
			allValues = new DoubleList[squares.length];
			for( i = 0; i < allValues.length; i++ ) {
				allValues[i] = new DoubleList();
			}
		}
		NumericalResult curr1, curr2;
		for (i = 0; i < results.length; i++) {
			curr1 = r1.getResultAt(i);
			curr2 = r2.getResultAt(i);
			// check name, ...
			if (!curr1.isComparableResult(curr2)) {
				throw new AdditionImpossibleException();
			}
			squares[i] = r1.squares[i] + r2.squares[i];
			results[i] = new NumericalResult(curr1.getName(), curr1
					.getComment(), ((Double) curr1.getValue())
					+ ((Double) curr2.getValue()));
			if( aggregateAll ) {
				allValues[i].addAll( r1.allValues[i] );
				allValues[i].addAll( r2.allValues[i] );
			}
		}
		MeanResultSet mrs = new MeanResultSet(results, r1.infos, squares, r1.count
				+ r2.count, aggregateAll, allValues );
		return mrs;
	}

	/**
	 * Adds {@link NumericalResultSet}s to this {@link MeanResultSet}. The
	 * {@link NumericalResultSet}s are handled as one result. So if you call
	 * this method with e.g. three arguments this is the same as adding one
	 * combing argument and not the same as calling the method three times, each
	 * time with one argument.
	 * 
	 * @param rs
	 *            the {@link NumericalResultSet}s to be added
	 * 
	 * @throws InconsistentResultNumberException
	 *             if the number of results differ
	 * @throws IllegalValueException
	 *             if the new (merged) value could not be set
	 * @throws AdditionImpossibleException
	 *             if some results are not comparable (name, comment, type)
	 */
	public synchronized void addResults(NumericalResultSet... rs)
			throws InconsistentResultNumberException, IllegalValueException,
			AdditionImpossibleException {
		int anz = 0, i = 0, idx = 0;
		for (; i < rs.length; i++) {
			if (rs[i] != null) {
				anz += rs[i].getNumberOfResults();
			}
		}

		boolean first = count == 0;
		if ( first ) {
			results = new AnnotatedEntityList<Result>( anz );
			squares = new double[anz];
			if( aggregateAll ) {
				allValues = new DoubleList[squares.length];
				for( i = 0; i < allValues.length; i++ ) {
					allValues[i] = new DoubleList();
				}
			}
		} else if (anz != this.getNumberOfResults()) {
			throw new InconsistentResultNumberException();
		}

		NumericalResult curr;
		for (int r = 0; r < rs.length; r++) {
			if (rs[r] != null) {
				int end = rs[r].getNumberOfResults();
				for ( i = 0; i < end; i++, idx++) {
					curr = rs[r].getResultAt(i);
					double currVal = 0;
					if (curr.getDatatype() == DataType.DOUBLE) {
						currVal = (Double) curr.getValue();
					} else {
						currVal = (Integer) curr.getValue();
					}
					if( allValues != null ) {
						allValues[idx].add( currVal );
					}
					squares[idx] += currVal * currVal;
					if ( first ) {
						results.add( new NumericalResult(curr.getName(), curr.getComment(), currVal) );
					} else {
						if (!results.get( idx ).isCastableResult(curr)) {
							throw new AdditionImpossibleException();
						} else {
							((NumericalResult) results.get( idx )).setResult( currVal + ((Number) results.get( idx ).getValue()).doubleValue() );
						}
					}
				}
			}
		}
		count++;
	}

	/**
	 * Returns the means and (if possible the) standard errors of the results in
	 * this {@link MeanResultSet} as a new {@link NumericalResultSet}.
	 * 
	 * @return the means and (if possible the) standard errors of the results in
	 *         this {@link MeanResultSet}
	 */
	public NumericalResultSet getStatistics() {
		int factor;
		if (count > 1) {
			factor = 2;
		} else {
			factor = 1;
		}
		NumericalResult[] resultsTemp = new NumericalResult[results.size()
				* factor];
		double n = count;
		for (int i = 0; i < results.size(); i++) {
			resultsTemp[i * factor] = new NumericalResult(results.get( i ).getName(),
					results.get( i ).getComment(), ((Double) results.get( i ).getValue())
							/ n);
			if (count > 1) {
				resultsTemp[(i * factor) + 1] = new NumericalResult(
						"Standard error of " + results.get( i ).getName(),
						"Standard error of the values of "
								+ results.get( i ).getName(), Math
								.sqrt((squares[i] / n - ((Double) resultsTemp[i
										* factor].getValue())
										* ((Double) resultsTemp[i * factor]
												.getValue()))
										/ (n - 1)));
			}
		}

		return new NumericalResultSet(resultsTemp);
	}
	
	/**
	 * Returns all values of the result with index <code>index</code> if available otherwise <code>null</code>.
	 * 
	 * @param index the index of the result
	 * 
	 * @return all values if available otherwise <code>null</code>
	 */
	public double[] getAllValues( int index ) {
		if( allValues == null ) {
			return null;
		} else {
			return allValues[index].toArray();
		}
	}

	private void setCount() {
		infos[infos.length - 1] = new NumericalResult(
				"evaluations",
				"the number of different results, each coming from one iteration of a crossvalidation",
				count);
	}

	/**
	 * Returns some information for this {@link MeanResultSet}.
	 * 
	 * @return the information for this {@link MeanResultSet}
	 */
	public ResultSet getInfos() {
		setCount();
		return new ResultSet(infos);
	}

	/**
	 * Class for the exception that is thrown if a {@link NumericalResultSet} is
	 * added to the {@link MeanResultSet} that has a number of results which is
	 * not equal to the number of results of the previously added results.
	 * 
	 * @author Jan Grau
	 */
	public static class InconsistentResultNumberException extends Exception {
		private static final long serialVersionUID = 1L;

		/**
		 * Constructs a new {@link InconsistentResultNumberException} with an
		 * appropriate error message.
		 * 
		 */
		public InconsistentResultNumberException() {
			super("Number of results differs.");
		}

	}

	/**
	 * Class for the exception that is thrown if two {@link MeanResultSet}s
	 * should be added that do not match.
	 * 
	 * @author Jens Keilwagen
	 */
	public static class AdditionImpossibleException extends Exception {
		private static final long serialVersionUID = 1L;

		/**
		 * Constructs a new {@link AdditionImpossibleException} with an
		 * appropriate error message.
		 * 
		 */
		public AdditionImpossibleException() {
			super("The addition is impossible, since the objects do not match.");
		}
	}
}
