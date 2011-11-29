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
import de.jstacs.classifier.scoringFunctionBased.gendismix.LearningPrinciple;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;

/**
 * This class implements a {@link CollectionParameter} based on an {@link Enum}.
 * Internally it is based on {@link String}s, i.e. the names of the {@link Enum}
 * constants. The methods {@link #setDefault(Object)} and
 * {@link #setValue(Object)} can be used with {@link String}s or with the
 * {@link Enum} constants.
 * 
 * @author Jens Keilwagen
 */
public class EnumParameter extends CollectionParameter {
	private Class<? extends Enum> enumInstance;
	private Enum[] enumConstants;

	private static String[] getKeys(Class<? extends Enum> e) {
		Enum[] all = e.getEnumConstants();
		String[] keys = new String[all.length];
		for (int i = 0; i < all.length; i++) {
			keys[i] = all[i].name();
		}
		return keys;
	}

	/**
	 * The main constructor.
	 * 
	 * @param enumInstance
	 *            the {@link Enum} class, e.g.
	 *            {@link de.jstacs.data.Sample.PartitionMethod}.class
	 * @param comment
	 *            a comment on this parameter
	 * @param required
	 *            <code>true</code> if this {@link EnumParameter} is required,
	 *            <code>false</code> otherwise
	 * 
	 * @throws ParameterException
	 *             is never thrown but exists due to the class hierarchy
	 */
	public EnumParameter(Class<? extends Enum> enumInstance, String comment,
			boolean required) throws ParameterException {
		super(DataType.STRING, getKeys(enumInstance), getKeys(enumInstance),
				enumInstance.getSimpleName(), comment, required);
		this.enumInstance = enumInstance;
		enumConstants = enumInstance.getEnumConstants();
	}

	/**
	 * This constructor creates an instance and set the default value.
	 * 
	 * @param enumInstance
	 *            the {@link Enum} class, e.g.
	 *            {@link de.jstacs.data.Sample.PartitionMethod}.class
	 * @param comment
	 *            a comment on this parameter
	 * @param required
	 *            <code>true</code> if this {@link EnumParameter} is required,
	 *            <code>false</code> otherwise
	 * @param defaultValue the default value of this parameter
	 * 
	 * @throws ParameterException
	 *             is never thrown but exists due to the class hierarchy
	 *
	 * @see EnumParameter#EnumParameter(Class, String, boolean)
	 * @see EnumParameter#setDefault(Object)
	 */
	public EnumParameter(Class<? extends Enum> enumInstance, String comment,
			boolean required, String defaultValue ) throws ParameterException {
		this(enumInstance,comment,required);
		setDefault( defaultValue );
	}

	
	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Restores an instance of {@link EnumParameter} from a XML representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link StringBuffer} <code>representation</code> could
	 *             not be parsed
	 */
	public EnumParameter(StringBuffer representation)
			throws NonParsableException {
		super(representation);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @seede.jstacs.parameters.CollectionParameter#appendCollection(java.lang.
	 * StringBuffer)
	 */
	@Override
	protected void appendCollection(StringBuffer buf) {
		XMLParser.appendObjectWithTags(buf, super.getValue().toString(),
				"selectedEnum");
		if (hasDefault()) {
			XMLParser.appendObjectWithTags(buf, parameters.getParameterAt(
					getDefault()).getValue().toString(), "defaultSelectedEnum");
		}
		XMLParser.appendObjectWithTags(buf, enumInstance.getName(), "enumName");
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.parameters.CollectionParameter#extractCollection(java.lang.
	 * StringBuffer)
	 */
	@Override
	protected void extractCollection(StringBuffer buf)
			throws NonParsableException {
		try {
			enumInstance = (Class<? extends Enum>) Class.forName(XMLParser
					.extractObjectForTags(buf, "enumName", String.class ));// TODO XMLP14CONV This and (possibly) the following lines have been converted automatically
			enumConstants = enumInstance.getEnumConstants();
			createParameterSet(getKeys(enumInstance), getKeys(enumInstance),
					null);
			if (hasDefault()) {
				setDefault(XMLParser.extractObjectForTags(buf, "defaultSelectedEnum", String.class ));// TODO XMLP14CONV This and (possibly) the following lines have been converted automatically
			}
			setValue(XMLParser.extractObjectForTags(buf, "selectedEnum", String.class ));
		} catch (Exception e) {
			throw new NonParsableException(e.getMessage());
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.CollectionParameter#getValue()
	 */
	@Override
	public Enum getValue() {
		// TODO return ...valueOf( (String) super.getValue() );
		int idx = 0;
		Object v = super.getValue();
		while (!enumConstants[idx].name().equals(v)) {
			idx++;
		}
		return enumConstants[idx];
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.CollectionParameter#setValue(java.lang.Object)
	 */
	@Override
	public void setValue(Object value) throws IllegalValueException {
		if (value != null && value instanceof Enum) {
			if (!enumInstance.isInstance(value)) {
				throw new IllegalValueException("Wrong Enum type.");
			} else {
				value = ((Enum) value).name();
			}
		}
		super.setValue(value);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.parameters.CollectionParameter#setDefault(java.lang.Object)
	 */
	@Override
	public void setDefault(Object defaultValue) throws IllegalValueException {
		if (defaultValue != null && defaultValue instanceof Enum) {
			if (!enumInstance.isInstance(defaultValue)) {
				throw new IllegalValueException("Wrong Enum type.");
			} else {
				defaultValue = ((Enum) defaultValue).name();
			}
		}
		super.setDefault(defaultValue);
	}
}
