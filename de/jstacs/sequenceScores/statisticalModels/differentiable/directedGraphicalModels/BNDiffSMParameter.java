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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Jstacs. If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package de.jstacs.sequenceScores.statisticalModels.differentiable.directedGraphicalModels;

import de.jstacs.Storable;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;

/**
 * Class for the parameters of a {@link BayesianNetworkDiffSM}. Each
 * parameter holds its current value, the symbol and position it is responsible
 * for and the context of the current symbol. A parameter can either be free or
 * determined by other parameters (on the same simplex).
 * 
 * @author Jan Grau
 */
public class BNDiffSMParameter implements Storable, Cloneable {

	/**
	 * The current value of this parameter.
	 */
	private double value;
	/**
	 * The symbol (out of some {@link de.jstacs.data.alphabets.Alphabet}) this parameter
	 * is responsible for.
	 */
	protected byte symbol;
	/**
	 * The position of <code>symbol</code> this parameter is responsible for.
	 */
	protected int position;
	/**
	 * The context of this parameter. The length <code>context.length</code> is
	 * equal to the number of parents, <code>context[i]</code> holds the
	 * information for parent no <code>i</code>, <code>context[i][0]</code>
	 * holds the position of parent <code>i</code>, <code>context[i][1]</code>
	 * ,... hold the possible configurations of parent <code>i</code> in the
	 * context. Normally, only <code>context[i][1]</code> exists, i.e. only one
	 * configuration is allowed in the context.
	 */
	protected int[][] context;
	/**
	 * The counts for this parameter. Used for determination of plug-in
	 * parameters.
	 */
	protected double count;
	/**
	 * The pseudocount for this parameter. Used for determination of plug-in
	 * parameters.
	 */
	protected double pseudoCount;
	private boolean free;

	private int index;
	private double expValue;
	private Double z;
	private Double t;

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#clone()
	 */
	@Override
	public BNDiffSMParameter clone() throws CloneNotSupportedException {
		BNDiffSMParameter clone = (BNDiffSMParameter) super.clone();
		for (int i = 0; i < context.length; i++) {
			clone.context[i] = context[i].clone();
		}
		return clone;
	}

	/**
	 * Creates a new {@link BNDiffSMParameter}, that is {@link BNDiffSMParameter} no
	 * <code>index</code> in the list of {@link BNDiffSMParameter}s of the
	 * {@link BayesianNetworkDiffSM} and responsible for
	 * <code>symbol</code> at position <code>position</code> and pseudo count
	 * <code>pseudoCount</code>.
	 * 
	 * @param index
	 *            the index in the list of all parameters
	 * @param symbol
	 *            the symbol this parameter is responsible for
	 * @param position
	 *            the position of the symbol
	 * @param pseudoCount
	 *            the pseudocount
	 * @param free
	 *            indicates if this parameter is a free parameter
	 */
	public BNDiffSMParameter(int index, byte symbol, int position, double pseudoCount,
			boolean free) {
		this(index, symbol, position, new int[0][2], pseudoCount, free);
	}

	/**
	 * Creates a new {@link BNDiffSMParameter}, that is {@link BNDiffSMParameter} no
	 * <code>index</code> in the list of {@link BNDiffSMParameter}s of the
	 * {@link BayesianNetworkDiffSM} and responsible for
	 * <code>symbol</code> at position <code>position</code> having context
	 * <code>context</code> and pseudocount <code>pseudoCount</code>.
	 * 
	 * @param index
	 *            the index in the list of all parameters
	 * @param symbol
	 *            the symbol this parameter is responsible for
	 * @param position
	 *            the position of the symbol
	 * @param context
	 *            the context of this parameter
	 * @param pseudoCount
	 *            the pseudocount
	 * @param free
	 *            indicates if this parameter is a free parameter
	 */
	public BNDiffSMParameter(int index, byte symbol, int position, int[][] context,
			double pseudoCount, boolean free) {
		this.index = index;
		this.symbol = symbol;
		this.position = position;
		this.context = context;
		this.pseudoCount = pseudoCount;
		this.count = pseudoCount;
		this.value = 0;
		this.expValue = 1;
		this.free = free;
		this.z = null;
		this.t = null;
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Re-creates a parameter from its XML-representation as returned by the
	 * method {@link #toXML()}.
	 * 
	 * @param representation
	 *            the XML code as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link StringBuffer} <code>representation</code> could
	 *             not be parsed
	 */
	public BNDiffSMParameter(StringBuffer representation) throws NonParsableException {
		representation = XMLParser.extractForTag(representation, "parameter");
		value = XMLParser.extractObjectForTags(representation, "value", double.class );
		expValue = Math.exp(value);
		index = XMLParser.extractObjectForTags(representation, "index", int.class );
		pseudoCount = XMLParser.extractObjectForTags(representation, "pseudoCount", double.class );
		symbol = XMLParser.extractObjectForTags(representation, "symbol", byte.class );
		position = XMLParser.extractObjectForTags(representation, "position", int.class );
		context = XMLParser.extractObjectForTags(representation, "context", int[][].class );
		count = XMLParser.extractObjectForTags(representation, "count", double.class );
		free = XMLParser.extractObjectForTags(representation, "free", boolean.class );
		z = XMLParser.extractObjectForTags(representation, "z", Double.class );
		t = XMLParser.extractObjectForTags(representation, "t", Double.class );
		/*
		String temp = XMLParser.extractObjectForTags(representation, "z", String.class );
		z = temp.equals("null") ? null : Double.parseDouble(temp);
		temp = XMLParser.extractObjectForTags(representation, "t", String.class );
		t = temp.equals("null") ? null : Double.parseDouble(temp);
		*/
	}

	/**
	 * Returns the pseudocount as given in the constructor.
	 * 
	 * @return the pseudocount
	 * 
	 * @see BNDiffSMParameter#BNDiffSMParameter(int, byte, int, double, boolean)
	 * @see BNDiffSMParameter#BNDiffSMParameter(int, byte, int, int[][], double, boolean)
	 */
	public double getPseudoCount() {
		return pseudoCount;
	}

	/**
	 * Resets the counts to the pseudocounts and the value to <code>0</code>.
	 */
	public void reset() {
		this.count = pseudoCount;
		this.value = 0;
		this.expValue = 1;
	}

	/**
	 * Returns the depth of the tree, i.e. the number of parents of this
	 * parameter.
	 * 
	 * @return the depth of the tree
	 */
	public int getDepth() {
		return context.length;
	}

	/**
	 * Prints the counts and the value of this parameter to {@link System#out}.
	 */
	public void print() {
		System.out.println(symbol + " c: " + count);
		System.out.println(symbol + ": " + value);
	}

	/**
	 * Indicates if the {@link Sequence} <code>seq</code> fulfills all
	 * requirements defined in the {@link #context}.
	 * 
	 * @param seq
	 *            the sequence
	 * 
	 * @return <code>true</code> if this parameter is responsible for
	 *         <code>seq</code>, <code>false</code> otherwise
	 */
	public double doesApplyFor(Sequence seq) {
		if (seq.discreteVal(position) == symbol) {
			for (int i = 0; i < context.length; i++) {
				boolean app = false;
				for (int j = 1; j < context[i].length; j++) {
					if (seq.discreteVal(context[i][0]) == context[i][j]) {
						app = true;
					}
				}
				if (!app) {
					return 0.0;
				}
			}
			return 1.0;
		} else {
			return 0.0;
		}
	}

	/**
	 * Returns the current value of this parameter.
	 * 
	 * @return the current value
	 */
	public double getValue() {
		return value;
	}

	/**
	 * Sets the current value of this parameter.
	 * 
	 * @param value
	 *            the new value to be set
	 */
	public void setValue(double value) {
		this.value = value;
		this.expValue = Math.exp(value);
	}

	/**
	 * Returns <code>Math.exp({@link #getValue()})</code>, which is pre-computed.
	 * 
	 * @return the exponential value of this parameter
	 */
	public double getExpValue() {
		return this.expValue;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.Storable#toXML()
	 */
	public StringBuffer toXML() {
		StringBuffer buf = new StringBuffer();
		XMLParser.appendObjectWithTags(buf, value, "value");
		XMLParser.appendObjectWithTags(buf, symbol, "symbol");
		XMLParser.appendObjectWithTags(buf, index, "index");
		XMLParser.appendObjectWithTags(buf, pseudoCount, "pseudoCount");
		XMLParser.appendObjectWithTags(buf, position, "position");
		XMLParser.appendObjectWithTags(buf, context, "context");
		XMLParser.appendObjectWithTags(buf, count, "count");
		XMLParser.appendObjectWithTags(buf, free, "free");
		
		XMLParser.appendObjectWithTags(buf, z, "z");
		XMLParser.appendObjectWithTags(buf, t, "t");
		/*
		if (z != null) {
			XMLParser.appendObjectWithTags(buf, z, "z");
		} else {
			XMLParser.appendObjectWithTags(buf, null, "z");
		}
		if (t != null) {
			XMLParser.appendObjectWithTags(buf, t, "t");
		} else {
			XMLParser.appendObjectWithTags(buf, "null", "t");
		}
		/**/
		XMLParser.addTags(buf, "parameter");
		return buf;
	}

	/**
	 * Resets all internal normalization constants to <code>null</code>.
	 */
	public void invalidateNormalizers() {
		this.z = null;
		this.t = null;
	}

	/**
	 * Returns the partial derivative of the normalization constant with respect
	 * to this parameter.
	 * 
	 * @return the partial derivative
	 * 
	 * @throws Exception
	 *             if no normalization constants have been pre-computed
	 */
	public double getLogPartialNormalizer() throws Exception {
		if (z == null || t == null) {
			throw new Exception("No valid normalizers available for parameter "
					+ index + " at position " + position + ": z=" + z + ", t="
					+ t + ".");
		}
		return t + value + z;
	}

	/**
	 * Returns the current counts for this parameter.
	 * 
	 * @return the current counts
	 */
	public double getCounts() {
		return count;
	}

	/**
	 * Adds <code>count2</code> to the counts of this parameter.
	 * 
	 * @param count2
	 *            the counts to be added
	 */
	public void addCount(double count2) {
		this.count += count2;
	}

	/**
	 * Indicates if this parameter is free.
	 * 
	 * @return <code>true</code> if the parameter is free, <code>false</code>
	 *         otherwise
	 */
	public boolean isFree() {
		return free;
	}

	/**
	 * Returns the position of the symbol this parameter is responsible for as
	 * defined in the constructor.
	 * 
	 * @return the position of this parameter this parameter is responsible for
	 * 
	 * @see BNDiffSMParameter#BNDiffSMParameter(int, byte, int, double, boolean)
	 * @see BNDiffSMParameter#BNDiffSMParameter(int, byte, int, int[][], double, boolean)
	 */
	public int getPosition() {
		return position;
	}

	/**
	 * Returns the index of this parameter as defined in the constructor.
	 * 
	 * @return the index of this parameter in the list of all parameters
	 * 
	 * @see BNDiffSMParameter#BNDiffSMParameter(int, byte, int, double, boolean)
	 * @see BNDiffSMParameter#BNDiffSMParameter(int, byte, int, int[][], double, boolean)
	 */
	public int getIndex() {
		return index;
	}

	/**
	 * Sets the part of the normalization constant of parameters before this
	 * parameter in the structure of the network.
	 * 
	 * @param t
	 *            the normalization constant
	 */
	void setLogT(Double t) {
		this.t = t;
	}

	/**
	 * Sets the part of the normalization constant of parameters after this
	 * parameter in the structure of the network.
	 * 
	 * @param z
	 *            the normalization constant
	 */
	void setLogZ(Double z) {
		this.z = z;
	}

	/**
	 * Returns the part of the normalization constant of parameters after this
	 * parameter in the structure of the network.
	 * 
	 * @return the normalization constant
	 */
	public double getLogZ() {
		return z;
	}

	/**
	 * Returns the part of the normalization constant of parameters before this
	 * parameter in the structure of the network.
	 * 
	 * @return the normalization constant
	 */
	public double getLogT() {
		return t;
	}

}
