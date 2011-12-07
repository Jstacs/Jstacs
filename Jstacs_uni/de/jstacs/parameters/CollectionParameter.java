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
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.SimpleParameter.DatatypeNotValidException;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;
import de.jstacs.utils.galaxy.GalaxyAdaptor;

/**
 * Class for a collection parameter, i.e. a parameter that provides some
 * collection of possible values the user can choose from. Instances of this
 * class can be used as parameters in a {@link ParameterSet}.
 * 
 * @see de.jstacs.parameters.ParameterSet
 * 
 * @author Jan Grau
 */
public class CollectionParameter extends Parameter implements Rangeable, GalaxyConvertible {

	/**
	 * The internal {@link ParameterSet} that holds the possible values
	 */
	protected ParameterSet parameters;

	/**
	 * The number of the currently selected value in <code>parameters</code>
	 */
	private int selected;

	/**
	 * The number of the option selected by default
	 */
	private int defaultSelected;

	/**
	 * <code>true</code> if the user has selected an item
	 */
	protected boolean userSelected;

	/**
	 * <code>true</code> if the current parameter is required, false otherwise
	 */
	private boolean required;

	/**
	 * If a value was illegal for the collection parameter, this field holds the
	 * error message.
	 */
	protected String errorMessage;

	/**
	 * Indicates if this {@link CollectionParameter} shall be rangeable
	 */
	private boolean rangeable;

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#clone()
	 */
	@Override
	public CollectionParameter clone() throws CloneNotSupportedException {
		CollectionParameter clone = (CollectionParameter) super.clone();
		clone.parameters = parameters == null ? null : parameters.clone();
		return clone;
	}

	// default constructor
	private CollectionParameter(DataType datatype, String name, String comment,
			boolean required) {
		this.datatype = datatype;
		this.comment = comment;
		this.name = name;
		this.required = required;
		this.userSelected = false;
		this.rangeable = true;
		this.defaultSelected = -1;
		this.selected = 0;
	}

	/**
	 * Creates a new {@link CollectionParameter} from the necessary field. This
	 * constructor should be used to clone a current instance.
	 * 
	 * @param options
	 *            the options of the {@link CollectionParameter}
	 * @param selected
	 *            the currently selected value
	 * @param defaultSelected
	 *            the value selected by default
	 * @param userSelected
	 *            <code>true</code> if the current value was selected by the
	 *            user, <code>false</code> otherwise
	 * @param name
	 *            the name of the parameter
	 * @param comment
	 *            a comment on the parameter
	 * @param required
	 *            <code>true</code> if this {@link CollectionParameter} is
	 *            required, <code>false</code> otherwise
	 * @param datatype
	 *            the data type of the parameters in the collection
	 * @param errorMessage
	 *            the error message of the last error or <code>null</code>
	 * @param rangeable
	 *            <code>true</code> if the current instance is rangeable
	 */
	protected CollectionParameter(ParameterSet options, int selected,
			int defaultSelected, boolean userSelected, String name,
			String comment, boolean required, DataType datatype,
			String errorMessage, boolean rangeable) {
		this.parameters = options;
		this.selected = selected;
		this.defaultSelected = defaultSelected;
		this.userSelected = userSelected;
		this.name = name;
		this.comment = comment;
		this.required = required;
		this.datatype = datatype;
		this.errorMessage = errorMessage;
		this.rangeable = rangeable;
	}

	/**
	 * Constructor for a {@link CollectionParameter}.
	 * 
	 * @param datatype
	 *            the data type of the parameters in the collection
	 * @param keys
	 *            the keys/names of the values in the collection, this is the
	 *            name the user will see in the user interface
	 * @param values
	 *            the values the names stand for, this array must be of the same
	 *            length as <code>keys</code>, a key at a certain position
	 *            belongs to the value at the same position in the array
	 * @param name
	 *            the name of the parameter
	 * @param comment
	 *            a comment on the parameter
	 * @param required
	 *            <code>true</code> if this {@link CollectionParameter} is
	 *            required, <code>false</code> otherwise
	 * 
	 * @throws InconsistentCollectionException
	 *             if the length of <code>keys</code> and the
	 *             <code>values</code> is different or the collection is
	 *             inconsistent for some other reason
	 * @throws IllegalValueException
	 *             if one of the values in <code>values</code> is not of type
	 *             <code>datatype</code>
	 * @throws DatatypeNotValidException
	 *             if the <code>datatype</code> is not one of the allowed values
	 * 
	 * @see CollectionParameter#CollectionParameter(DataType, String[],
	 *      Object[], String[], String, String, boolean)
	 */
	public CollectionParameter(DataType datatype, String[] keys,
			Object[] values, String name, String comment, boolean required)
			throws InconsistentCollectionException, IllegalValueException,
			DatatypeNotValidException {
		this(datatype, keys, values, null, name, comment, required);
	}

	/**
	 * Constructor for a {@link CollectionParameter}.
	 * 
	 * @param datatype
	 *            the data type of the parameters in the collection
	 * @param keys
	 *            the keys/names of the values in the collection, this is the
	 *            name the user will see in the user interface
	 * @param values
	 *            the values the names stand for, this array must be of the same
	 *            length as <code>keys</code>, a key at a certain position
	 *            belongs to the value at the same position in the array
	 * @param comments
	 *            the comments on the possible values
	 * @param name
	 *            the name of the parameter
	 * @param comment
	 *            a comment on the parameter
	 * @param required
	 *            <code>true</code> if the {@link CollectionParameter} is
	 *            required, <code>false</code> otherwise
	 * 
	 * @throws InconsistentCollectionException
	 *             if the length of <code>keys</code> and <code>values</code> is
	 *             different or the collection is inconsistent for some other
	 *             reason
	 * @throws IllegalValueException
	 *             if one of the values in <code>values</code> is not of type
	 *             <code>datatype</code>
	 * @throws DatatypeNotValidException
	 *             if the <code>datatype</code> is not one of the allowed values
	 * 
	 * @see CollectionParameter#createParameterSet(Object[], String[], String[])
	 */
	public CollectionParameter(DataType datatype, String[] keys,
			Object[] values, String[] comments, String name, String comment,
			boolean required) throws InconsistentCollectionException,
			IllegalValueException, DatatypeNotValidException {
		this(datatype, name, comment, required);

		if (keys.length != values.length
				|| (comments != null && keys.length != comments.length)) {
			throw new InconsistentCollectionException(
					"You have to define the same number of keys and values for a CollectionParameter!");
		}

		createParameterSet(values, keys, comments);
	}

	/**
	 * Constructor for a {@link CollectionParameter} from an array of
	 * {@link ParameterSet}s. This constructor can be used to easily construct a
	 * {@link CollectionParameter} that lets the user select from a list of
	 * possible options that all require an own set of {@link Parameter}s to be
	 * instantiated.
	 * 
	 * @param values
	 *            the array of {@link ParameterSet}s
	 * @param keys
	 *            the keys/names of the values in the collection, this is the
	 *            name the user will see in the user interface
	 * @param comments
	 *            the comments on the possible values
	 * @param name
	 *            the name of the parameter
	 * @param comment
	 *            a comment on the parameter
	 * @param required
	 *            <code>true</code> if the parameter is required,
	 *            <code>false</code> otherwise
	 * 
	 * @see CollectionParameter#createParameterSet(Object[], String[], String[])
	 */
	public CollectionParameter(ParameterSet[] values, String[] keys,
			String[] comments, String name, String comment, boolean required) {
		this(DataType.PARAMETERSET, name, comment, required);
		try {
			createParameterSet(values, keys, comments);
		} catch (Exception doesnothappen) {
			doesnothappen.printStackTrace();
		}
	}

	/**
	 * Constructor for a {@link CollectionParameter} from an array of
	 * {@link ParameterSet}s. This constructor can be used to easily construct a
	 * {@link CollectionParameter} that lets the user select from a list of
	 * possible options that all require an own set of {@link Parameter}s to be
	 * instantiated.
	 * 
	 * @param values
	 *            the array of {@link ParameterSet}s
	 * @param name
	 *            the name of the parameter
	 * @param comment
	 *            a comment on the parameter
	 * @param required
	 *            <code>true</code> if the parameter is required,
	 *            <code>false</code> otherwise
	 * 
	 * @see CollectionParameter#createParameterSet(Object[], String[], String[])
	 */
	public CollectionParameter(InstanceParameterSet[] values, String name,
			String comment, boolean required) {
		this(DataType.PARAMETERSET, name, comment, required);
		try {
			createParameterSet(values, null, null);
		} catch (Exception doesnothappen) {
		}
	}

	/**
	 * Creates a new {@link ParameterSet} from an array of values, an array of
	 * names and an array of comments.
	 * 
	 * @param values
	 *            the values the names stand for, this array must be of the same
	 *            length as <code>keys</code>, a key at a certain position
	 *            belongs to the value at the same position in the array
	 * @param keys
	 *            the keys/names of the values in the collection, this is the
	 *            name the user will see in the user interface
	 * @param comments
	 *            the comments on the possible values
	 * 
	 * @throws IllegalValueException
	 *             if one of the values in <code>values</code> is not of type
	 *             <code>datatype</code>
	 * @throws DatatypeNotValidException
	 *             if the data type is not allowed, i.e. is not one of the
	 *             primitive data types, {@link String} or
	 *             {@link DataType#PARAMETERSET}
	 */
	protected void createParameterSet(Object[] values, String[] keys,
			String[] comments) throws DatatypeNotValidException,
			IllegalValueException {
		Parameter[] pars = new Parameter[values.length];
		InstanceParameterSet p;
		String c, k;
		for (int i = 0; i < pars.length; i++) {

			if (values[i] instanceof InstanceParameterSet) {
				p = (InstanceParameterSet) values[i];
				if (keys == null || keys[i] == null) {
					k = p.getInstanceName();
				} else {
					k = keys[i];
				}
				if (comments == null) {
					c = p.getInstanceComment();
				} else {
					c = comments[i];
				}
				pars[i] = new ParameterSetContainer(k, c, p);
			} else {
				if (keys == null || keys[i] == null) {
					throw new IllegalArgumentException(
							"You have to state the key for entity " + i);
				}
				if (comments != null) {
					c = comments[i];
				} else {
					c = "";
				}
				if (values[i] instanceof ParameterSet) {
					pars[i] = new ParameterSetContainer(keys[i], c,
							(ParameterSet) values[i]);
				} else {
					pars[i] = new SimpleParameter(datatype, keys[i], c, false);
					pars[i].setValue(values[i]);
				}
			}
		}

		parameters = new SimpleParameterSet(pars);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#isAtomic()
	 */
	@Override
	public boolean isAtomic() {
		// return !(getValue() instanceof ParameterSet);
		for (int i = 0; i < parameters.getNumberOfParameters(); i++) {
			if (parameters.getParameterAt(i).getValue() instanceof ParameterSet) {
				return false;
			}
		}
		return true;
	}

	/**
	 * Returns <code>true</code>, if this {@link CollectionParameter} has a
	 * default value.
	 * 
	 * @return if this {@link CollectionParameter} has a default value
	 * 
	 * @see CollectionParameter#CollectionParameter(ParameterSet, int, int,
	 *      boolean, String, String, boolean, DataType, String, boolean)
	 * @see CollectionParameter#setDefault(Object)
	 */
	public boolean hasDefault() {
		return defaultSelected > -1;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#hasDefaultOrIsSet()
	 */
	@Override
	public boolean hasDefaultOrIsSet() {
		if (getValue() instanceof ParameterSet) {
			if (!((ParameterSet) getValue()).hasDefaultOrIsSet()) {
				if (((ParameterSet) getValue()).getErrorMessage() != null) {
					errorMessage = "Selected value has the following error: "
							+ ((ParameterSet) getValue()).getErrorMessage();
				}
				return false;
			} else {
				return true;
			}
		} else if (isSet() || selected == defaultSelected) {
			return true;
		} else {
			return false;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#isSet()
	 */
	@Override
	public boolean isSet() {
		return isUserSelected();
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Restores an instance of {@link CollectionParameter} from a XML
	 * representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link StringBuffer} <code>representation</code> could
	 *             not be parsed
	 * 
	 * @see CollectionParameter#fromXML(StringBuffer)
	 */
	public CollectionParameter(StringBuffer representation)
			throws NonParsableException {
		fromXML(representation);
	}

	/**
	 * Returns the possible values in this collection.
	 * 
	 * @return the values/parameters in this collection
	 */
	public ParameterSet getParametersInCollection() {
		return parameters;
	}

	/**
	 * Sets the value returned by {@link #isRangeable()} to
	 * <code>rangeable</code>.
	 * 
	 * @param rangeable
	 *            the new value
	 */
	public void setRangeable(boolean rangeable) {
		this.rangeable = rangeable;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Rangeable#isRangeable()
	 */
	public boolean isRangeable() {
		return rangeable;
	}

	/**
	 * Returns <code>true</code> if the key specified by <code>value</code> is
	 * in the set of keys of this {@link CollectionParameter}.
	 * 
	 * @return <code>true</code> if the key specified by <code>value</code> is
	 *         in the set of keys of this {@link CollectionParameter},
	 *         <code>false</code> otherwise
	 */
	@Override
	public boolean checkValue(Object value) {
		if (value == null) {
			return false;
		}
		if (!(value instanceof String)) {
			errorMessage = "The value is not in the set of defined values (and not even a String)."
					+ value;
			return false;
		} else {
			for (int i = 0; i < parameters.getNumberOfParameters(); i++) {
				if (parameters.getParameterAt(i).getName().equals(value)) {
					errorMessage = null;
					return true;
				}
			}
		}
		errorMessage = "The value is not in the set of defined values: "
				+ value + ".";
		return false;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#getErrorMessage()
	 */
	@Override
	public String getErrorMessage() {
		return errorMessage;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#toXML()
	 */
	@Override
	public StringBuffer toXML() {
		StringBuffer buf = super.toXML();
		XMLParser.addTags(buf, "superParameter");
		XMLParser.appendObjectWithTags(buf, datatype, "datatype");
		XMLParser.appendObjectWithTags(buf, name, "name");
		XMLParser.appendObjectWithTags(buf, comment, "comment");
		XMLParser.appendObjectWithTags(buf, required, "required");
		XMLParser.appendObjectWithTags(buf, userSelected, "userSelected");
		XMLParser.appendObjectWithTags(buf, errorMessage, "errorMessage");
		XMLParser.appendObjectWithTags(buf, selected, "selected");
		XMLParser.appendObjectWithTags(buf, defaultSelected, "defaultSelected");
		XMLParser.appendObjectWithTags(buf, rangeable, "rangeable");
		appendCollection(buf);

		XMLParser.addTags(buf, "collectionParameter");

		return buf;
	}

	/**
	 * Appends the internal {@link ParameterSet} in its XML representation (
	 * {@link ParameterSet#toXML()}) to the {@link StringBuffer}
	 * <code>buf</code>.
	 * 
	 * @param buf
	 *            the {@link StringBuffer} this method appends to
	 */
	protected void appendCollection(StringBuffer buf) {
		XMLParser.appendObjectWithTags(buf, parameters, "collection");
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#fromXML(java.lang.StringBuffer)
	 */
	@Override
	protected void fromXML(StringBuffer representation)
			throws NonParsableException {
		representation = XMLParser.extractForTag(representation,"collectionParameter");
		super.fromXML(XMLParser.extractForTag(representation,"superParameter"));
		datatype = XMLParser.extractObjectForTags(representation, "datatype", DataType.class );// TODO XMLP14CONV This and (possibly) the following lines have been converted automatically
		name = XMLParser.extractObjectForTags(representation, "name", String.class );
		comment = XMLParser.extractObjectForTags(representation, "comment", String.class );
		required = XMLParser.extractObjectForTags(representation, "required", boolean.class );
		userSelected = XMLParser.extractObjectForTags(representation, "userSelected", boolean.class );// TODO XMLP14CONV This and (possibly) the following lines have been converted automatically
		errorMessage = XMLParser.parseString( XMLParser.extractObjectForTags(representation, "errorMessage", String.class ) );
		selected = XMLParser.extractObjectForTags(representation, "selected", int.class );
		defaultSelected = XMLParser.extractObjectForTags(representation, "defaultSelected", int.class );
		StringBuffer help = XMLParser.extractForTag(representation, "rangeable");
		if (help == null) {
			rangeable = false;
		} else {
			rangeable = Boolean.parseBoolean(help.toString());
		}
		extractCollection(representation);
	}

	/**
	 * Extracts the internal {@link ParameterSet} from its XML representation (
	 * {@link ParameterSet#toXML()}). Reverse method to
	 * {@link CollectionParameter#appendCollection(StringBuffer)}.
	 * 
	 * @param buf
	 *            the {@link StringBuffer} containing the XML representation
	 * 
	 * @throws NonParsableException
	 *             if the XML code could not be parsed
	 */
	protected void extractCollection(StringBuffer buf)
			throws NonParsableException {
		parameters = new SimpleParameterSet(XMLParser.extractForTag(buf,"collection"));
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#isRequired()
	 */
	@Override
	public boolean isRequired() {
		return required;
	}

	/**
	 * Returns <code>true</code> if the option at position <code>idx</code> is
	 * selected.
	 * 
	 * @param idx
	 *            the position
	 * 
	 * @return <code>true</code> if the option at position <code>idx</code> is
	 *         selected, <code>false</code> otherwise
	 */
	public boolean isSelected(int idx) {
		return idx == selected;
	}
	
	/**
	 * Returns the index of the selected value.
	 * 
	 * @return the index of the selected value
	 */
	public int getSelected() {
		return selected;
	}

	/**
	 * Returns the index of the default selected value.
	 * 
	 * @return the index of the default selected value
	 */
	protected int getDefault() {
		return defaultSelected;
	}

	/**
	 * Returns <code>true</code> if the value was selected by the user.
	 * 
	 * @return <code>true</code> if the value was selected by the user.
	 */
	public boolean isUserSelected() {
		return userSelected;
	}

	/**
	 * Sets the selected value to the one that is specified by the key
	 * <code>value</code>.
	 * 
	 * @param value
	 *            the key of the desired value
	 */
	@Override
	public void setValue(Object value) throws IllegalValueException {
		if (value == null) {
			return;
		}
		Object val2 = value instanceof InstanceParameterSet ? ((InstanceParameterSet)value).getInstanceName() : value;
		if (checkValue(val2)) {
			for (int i = 0; i < parameters.getNumberOfParameters(); i++) {
				if (parameters.getParameterAt(i).getName().equals(val2)) {
					selected = i;
					userSelected = true;
					break;
				}
			}
			if( value != val2 ) {
				parameters.getParameterAt( selected ).setValue( value );
			}
		} else {
			String s = "", sep= "";
			for (int i = 0; i < parameters.getNumberOfParameters(); i++) {
				s += sep + parameters.getParameterAt(i).getName();
				sep=", ";
			}
			throw new IllegalValueException("Value (" + val2 + ") not in Collection (" + s + ")!");
		}

	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#setDefault(java.lang.Object)
	 */
	@Override
	public void setDefault(Object defaultValue) throws IllegalValueException {
		setValue(defaultValue);
		defaultSelected = selected;
		userSelected = false;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#simplify()
	 */
	@Override
	public void simplify() {
		for (int i = 0; i < parameters.getNumberOfParameters(); i++) {
			if (i != selected
					&& parameters.getParameterAt(i).getValue() instanceof ParameterSet) {
				parameters.getParameterAt(i).reset();
			} else {
				parameters.getParameterAt(i).simplify();
			}
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#reset()
	 */
	@Override
	public void reset() {
		selected = defaultSelected;
		simplify();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#getValue()
	 */
	@Override
	public Object getValue() {
		if (selected < parameters.getNumberOfParameters()) {
			return parameters.getParameterAt(selected).getValue();
		} else {
			return null;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#equals(java.lang.Object)
	 */
	@Override
	public boolean equals(Object o2) {
		if (o2 instanceof CollectionParameter) {
			ParameterSet parSet2 = ((CollectionParameter) o2)
					.getParametersInCollection();
			if (parSet2.getNumberOfParameters() != parameters
					.getNumberOfParameters()) {
				return false;
			} else {
				if (!(((CollectionParameter) o2).getName().equals(name) && ((CollectionParameter) o2)
						.getComment().equals(comment))) {
					return false;
				}

				for (int i = 0; i < parameters.getNumberOfParameters(); i++) {
					if (!parameters.getParameterAt(i).equals(
							parSet2.getParameterAt(i))) {
						return false;
					}
				}
				return true;
			}
		} else {
			return false;
		}
	}

	/**
	 * This exception is thrown if the {@link CollectionParameter} is
	 * inconsistent for some reason.
	 * 
	 * @author Jan Grau
	 */
	public class InconsistentCollectionException extends ParameterException {

		private static final long serialVersionUID = -2703514434545861722L;

		/**
		 * Constructs a new {@link InconsistentCollectionException} with message
		 * <code>message</code>.
		 * 
		 * @param message
		 *            the error message
		 */
		public InconsistentCollectionException(String message) {
			super(message);
		}

	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Rangeable#getRangedInstance()
	 */
	public Parameter getRangedInstance() throws Exception {
		boolean[] selected = new boolean[parameters.getNumberOfParameters()];
		selected[0] = true;
		MultiSelectionCollectionParameter par = new MultiSelectionCollectionParameter(
				this.parameters.clone(), selected, new boolean[parameters
						.getNumberOfParameters()], false, getName(),
				getComment(), required, datatype, errorMessage, 0, true);
		return par;
	}
	
	/*
	 * (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	public String toString() {
		String s = parameters.getParameterAt( 0 ).getName();//.getValue().toString();
		for( int i = 1; i < parameters.getNumberOfParameters(); i++ ) {
			s += ", " +  parameters.getParameterAt( i ).getName();//.getValue().toString();
		}
		return name + " (" + comment
			+ ", range={" + s + "}" 
			+ (defaultSelected>=0?", default = " + parameters.getParameterAt(defaultSelected).getValue():"")
			+ (required ? "" : ", OPTIONAL" )
			+ ")\t= " + getValue();
	}

	@Override
	public void toGalaxy( String namePrefix, String configPrefix, int depth, StringBuffer descBuffer, StringBuffer configBuffer ) throws Exception {
		StringBuffer buf = new StringBuffer();
		namePrefix = namePrefix+"_"+GalaxyAdaptor.getLegalName( getName() );
		
		for(int i=0;i<parameters.getNumberOfParameters();i++){
			if(isSelected( i )){
				XMLParser.appendObjectWithTagsAndAttributes( buf, parameters.getParameterAt( i ).getName() , "option", "value=\""+parameters.getParameterAt( i ).getName()+"\" selected=\"true\"", false );
			}else{
				XMLParser.appendObjectWithTagsAndAttributes( buf, parameters.getParameterAt( i ).getName() , "option", "value=\""+parameters.getParameterAt( i ).getName()+"\"", false );
			}
		}
		XMLParser.addTagsAndAttributes( buf, "param", "type=\"select\" format=\"text\" name=\""+namePrefix+"\" label=\""+getName()+"\" optional=\""+(!isRequired())+"\" help=\""+getComment()+"\"" );
	
		StringBuffer buf3 = new StringBuffer("${"+configPrefix+(isAtomic() ? "" : namePrefix+"_cond.")+namePrefix+"}");
		XMLParser.addTags( buf3, namePrefix );
		configBuffer.append( buf3 );
		if(!isAtomic()){
			for(int i=0;i<parameters.getNumberOfParameters();i++){
				StringBuffer temp = new StringBuffer();
				StringBuffer temp2 = new StringBuffer();
				((GalaxyConvertible)parameters.getParameterAt( i )).toGalaxy( namePrefix+"_opt"+i, configPrefix+namePrefix+"_cond.", depth+1, temp, temp2 );
				XMLParser.addTagsAndAttributes( temp, "when", "value=\""+parameters.getParameterAt( i ).getName()+"\"" );
				buf.append( temp );
				configBuffer.append( (i== 0 ? "#if " : "#elif ")+"$"+configPrefix+namePrefix+"_cond."+namePrefix+" == \""+parameters.getParameterAt( i ).getName()+"\"\n" );
				configBuffer.append( temp2 );
				configBuffer.append( "\n" );
			}
			configBuffer.append( "#end if\n" );
		}
		
		
		if(!isAtomic()){
			XMLParser.addTagsAndAttributes( buf, "conditional", "name=\""+namePrefix+"_cond\"" );
		}
		
		descBuffer.append( buf );
	}

	@Override
	public void fromGalaxy( String namePrefix, StringBuffer command ) throws Exception {
		namePrefix = namePrefix+"_"+GalaxyAdaptor.getLegalName( getName() );
		
		String selected = XMLParser.extractForTag( command, namePrefix ).toString();
		this.setValue( selected );
		if(this.getValue() instanceof GalaxyConvertible){
			((GalaxyConvertible)this.getValue()).fromGalaxy( namePrefix+"_opt"+getSelected(), command );
		}
	}
	
	
	
	
}