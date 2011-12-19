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

package de.jstacs.parameters;

import de.jstacs.DataType;
import de.jstacs.NonParsableException;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.AlphabetContainerParameterSet;
import de.jstacs.data.AlphabetContainer.AlphabetContainerType;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.validation.NumberValidator;

/**
 * Abstract class for a {@link ParameterSet} containing all parameters necessary
 * to construct an {@link Object} that implements
 * {@link de.jstacs.InstantiableFromParameterSet}. This parameter set handles the
 * {@link AlphabetContainer} and if necessary the length of a sequence, so it is
 * well suited as parameter set for {@link de.jstacs.models.AbstractModel} and
 * {@link de.jstacs.classifier.AbstractClassifier}.
 * 
 * @author Jan Grau, Jens Keilwagen
 * 
 * @see de.jstacs.InstantiableFromParameterSet
 * @see ParameterSet
 * @see de.jstacs.models.AbstractModel
 * @see de.jstacs.classifier.AbstractClassifier
 */
public abstract class SequenceScoringParameterSet extends InstanceParameterSet {

	/**
	 * The alphabet the model works on
	 */
	protected Parameter alphabet;

	/**
	 * The length of sequences the model can work on or <code>0</code> for
	 * arbitrary length
	 */
	protected Parameter length;

	/**
	 * <code>true</code> if the model can handle sequences of variable length
	 */
	private boolean variableLength;

	/**
	 * Constructs an {@link InstanceParameterSet} having empty parameter values.
	 * This constructor should only be used if the object can handle sequences
	 * of fixed length.
	 * 
	 * @param instanceClass
	 *            the class of the instance
	 * @param type
	 *            the type of the alphabet(s)
	 * @param simple
	 *            determines whether the alphabet should be simple
	 * 
	 * @see AlphabetContainerType
	 */
	public SequenceScoringParameterSet(Class instanceClass,
			AlphabetContainerType type, boolean simple) {
		this(instanceClass, type, simple, false);
	}

	/**
	 * Constructs a {@link SequenceScoringParameterSet} having empty parameter
	 * values. The user can specify a-priori if the object can handle sequences
	 * of variable lengths. If that is the case the object is not queried from
	 * the user as it is <code>0</code> anyway.
	 * 
	 * @param instanceClass
	 *            the (sub-)class of the instance
	 * @param type
	 *            the type of the alphabet(s)
	 * @param simple
	 *            determines whether the alphabet should be simple
	 * @param variableLength
	 *            <code>true</code> if the object can handle sequences of
	 *            arbitrary length, <code>false</code> otherwise
	 * 
	 * @see AlphabetContainerType
	 */
	public SequenceScoringParameterSet(Class instanceClass,
			AlphabetContainerType type, boolean simple, boolean variableLength) {
		super(instanceClass);
		this.variableLength = variableLength;
		try {
			this.alphabet = new ParameterSetContainer("Alphabet",
					"The alphabet the model works on",
					new AlphabetContainerParameterSet(type, simple));
			if (variableLength) {
				this.length = new SimpleParameter(DataType.INT, "Length",
						"The length of sequences the model can work on", true,
						new NumberValidator<Integer>(0, 0), 0);
			} else {
				this.length = new SimpleParameter(DataType.INT, "Length",
						"The length of sequences the model can work on", true,
						new NumberValidator<Integer>(1, Integer.MAX_VALUE));
			}
		} catch (Exception hopefullyDoesNotHappen) {
			RuntimeException ex = new RuntimeException( hopefullyDoesNotHappen );
			ex.setStackTrace( hopefullyDoesNotHappen.getStackTrace() );
			throw ex;
		}
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link SequenceScoringParameterSet} from its XML
	 * representation. Automatically calls the current implementation of
	 * {@link #fromXML(StringBuffer)}.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link StringBuffer} <code>representation</code> could
	 *             not be parsed
	 */
	public SequenceScoringParameterSet(StringBuffer representation)
			throws NonParsableException {
		super(representation);
	}

	/**
	 * Constructs a {@link SequenceScoringParameterSet} from an
	 * {@link AlphabetContainer} and the length of a sequence. This constructor
	 * can be used to implement a {@link SequenceScoringParameterSet} that is
	 * already instantiated with known parameter values.
	 * 
	 * @param instanceClass
	 *            the class of the instance
	 * @param alphabet
	 *            the {@link AlphabetContainer} for a sequence
	 * @param length
	 *            the length of sequence
	 * @param variableLength
	 *            <code>true</code> if the object can handle sequences of
	 *            arbitrary length, <code>false</code> otherwise
	 * 
	 * @throws Exception
	 *             if the alphabets or the length are not in the expected range
	 *             of values
	 */
	public SequenceScoringParameterSet(Class instanceClass,
			AlphabetContainer alphabet, int length, boolean variableLength)
			throws Exception {
		this(instanceClass, alphabet.getType(), alphabet.isSimple(),
				variableLength);
		this.length.setValue(length);
		this.alphabet.setValue(alphabet.getCurrentParameterSet());
	}

	/**
	 * Constructs a {@link SequenceScoringParameterSet} for an object that can
	 * handle sequences of variable length and with the
	 * {@link AlphabetContainer} <code>alphabet</code>. This constructor can be
	 * used to implement a {@link SequenceScoringParameterSet} that is already
	 * instantiated with known parameter values.
	 * 
	 * <br>
	 * 
	 * The length parameter is set to <code>0</code>.
	 * 
	 * @param instanceClass
	 *            the (sub-)class of the instance
	 * @param alphabet
	 *            the {@link AlphabetContainer} for a sequence
	 * 
	 * @throws Exception
	 *             if the alphabets or the length are not in the expected range
	 *             of values
	 * 
	 */
	public SequenceScoringParameterSet(Class instanceClass,
			AlphabetContainer alphabet) throws Exception {
		this(instanceClass, alphabet, 0, true);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.ParameterSet#hasDefaultOrIsSet()
	 */
	@Override
	public boolean hasDefaultOrIsSet() {
		boolean erg;
		erg = super.hasDefaultOrIsSet() && alphabet.hasDefaultOrIsSet()
				&& length.hasDefaultOrIsSet();

		if (erg) {
			// AlphabetContainer abc = getAlphabet();
			int l = 0;
			if (alphabet != null) {
				l = ((AlphabetContainerParameterSet) alphabet.getValue())
						.getPossibleLength();
			}
			erg &= l == 0 || ((Integer) length.getValue()).intValue() == l;
			if (!erg) {
				errorMessage = "The length of the alphabet and the length of the model must match!";
			}
		}
		return erg;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.ParameterSet#reset()
	 */
	@Override
	public void reset() {
		super.reset();
		length.reset();
		alphabet.reset();
	}

	/**
	 * Returns the {@link AlphabetContainer} of the current instance.
	 * 
	 * @return the {@link AlphabetContainer}
	 */
	public AlphabetContainer getAlphabetContainer() {
		try {
			AlphabetContainer cont = new AlphabetContainer(
					(AlphabetContainerParameterSet) alphabet.getValue());
			return cont;
		} catch (Exception e) {
			e.printStackTrace();
			return null;
		}
	}

	/**
	 * Returns the length of the sequences the model can work on.
	 * 
	 * @return the length of the sequences the model can work on
	 * 
	 * @throws IllegalArgumentException
	 *             if the length is not correct, i.e. the length is not suitable
	 *             for the chosen requirements
	 */
	public int getLength() throws IllegalArgumentException {
		int l = (Integer) length.getValue();
		if (variableLength && l != 0) {
			throw new IllegalArgumentException(
					"The model can handle sequences of variable length, but length is defined as "
							+ l);
		}
		return l;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.parameters.InstanceParameterSet#fromXML(java.lang.StringBuffer)
	 */
	@Override
	protected void fromXML(StringBuffer representation)
			throws NonParsableException {
		representation = XMLParser.extractForTag( representation, "sequenceScoringParameterSet" );
		super.fromXML(XMLParser.extractForTag( representation, "superParameters" ));
		variableLength = XMLParser.extractObjectForTags( representation, "variableLength", boolean.class );
		StringBuffer alphStringB = XMLParser.extractForTag( representation, "alphabet" );
		if (alphStringB == null) {
			alphabet = null;
		} else {
			alphabet = new ParameterSetContainer(alphStringB);
		}
		length = new SimpleParameter(XMLParser.extractForTag(representation,"length"));
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.ParameterSet#getNumberOfParameters()
	 */
	@Override
	public int getNumberOfParameters() {
		if (variableLength) {
			// System.out.println("homogeneous:
			// "+(super.getNumberOfParameters() + 1));
			return super.getNumberOfParameters() + 1;
		} else {
			// System.out.println("inhomogeneous:
			// "+(super.getNumberOfParameters() + 2));
			return super.getNumberOfParameters() + 2;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.ParameterSet#getParameterAt(int)
	 */
	@Override
	public Parameter getParameterAt(int i) {
		if (i < super.getNumberOfParameters()) {
			return super.getParameterAt(i);
		} else if (i == super.getNumberOfParameters()) {
			return alphabet;
		} else if (!variableLength) {
			return length;
		} else {
			throw new IndexOutOfBoundsException();
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.InstanceParameterSet#toXML()
	 */
	@Override
	public StringBuffer toXML() {
		StringBuffer buf = super.toXML();
		XMLParser.addTags( buf, "superParameters" );
		XMLParser.appendObjectWithTags(buf, variableLength, "variableLength");
		if (alphabet != null) {
			XMLParser.appendObjectWithTags(buf, alphabet, "alphabet");
		}
		XMLParser.appendObjectWithTags(buf, length, "length");
		XMLParser.addTags(buf, "sequenceScoringParameterSet");

		return buf;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#equals(java.lang.Object)
	 */
	@Override
	public boolean equals(Object o) {
		if (o instanceof SequenceScoringParameterSet) {
			SequenceScoringParameterSet comp = (SequenceScoringParameterSet) o;
			boolean erg = this.getInstanceClass().equals(
					comp.getInstanceClass())
					&& (getLength() == comp.getLength())
					&& (parameters.size() == comp.parameters.size())
					&& getAlphabetContainer().checkConsistency(comp.getAlphabetContainer());
			int i = 0;
			while (i < parameters.size() && erg) {
				erg &= parameters.get(i).equals(comp.parameters.get(i++));
			}
			return erg;
		} else {
			return false;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.ParameterSet#clone()
	 */
	@Override
	public SequenceScoringParameterSet clone()
			throws CloneNotSupportedException {
		SequenceScoringParameterSet res = (SequenceScoringParameterSet) super
				.clone();
		res.alphabet = alphabet.clone();
		res.length = length.clone();
		return res;
	}
}
