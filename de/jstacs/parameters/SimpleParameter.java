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
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.parameters.validation.ParameterValidator;
import de.jstacs.tools.ui.galaxy.GalaxyAdaptor;

/**
 * Class for a &quot;simple&quot; parameter. Simple parameters are those of the
 * primitive types and single {@link String}s.
 * 
 * @author Jan Grau
 */
public class SimpleParameter extends Parameter implements GalaxyConvertible {

	/**
	 * Determines if the parameter is required
	 */
	protected boolean required;

	/**
	 * The current value of the parameter
	 */
	protected Object value;

	/**
	 * The default value of the parameter
	 */
	protected Object defaultValue;

	/**
	 * The validator for the parameter values
	 */
	protected ParameterValidator validator;

	/**
	 * <code>true</code> if the parameters are set, <code>false</code> otherwise
	 */
	private boolean isSet;

	/**
	 * The error message or <code>null</code> if no error occurred
	 */
	private String errorMessage;

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs a {@link SimpleParameter} out of an XML representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link SimpleParameter} could not be restored from the
	 *             {@link StringBuffer} <code>representation</code>
	 */
	public SimpleParameter( StringBuffer representation ) throws NonParsableException {
		super(representation);
	}

	/**
	 * Constructor for a {@link SimpleParameter} without
	 * {@link ParameterValidator}.
	 * 
	 * @param datatype
	 *            the data type of the parameter value
	 * @param name
	 *            the name of the parameter
	 * @param comment
	 *            a comment on the parameter that tells the user some details
	 *            about it
	 * @param required
	 *            determines if the parameter is required
	 * 
	 * @throws DatatypeNotValidException
	 *             if <code>datatype</code> is not one of the allowed
	 *             {@link DataType}s
	 */
	public SimpleParameter(DataType datatype, String name, String comment,
			boolean required) throws DatatypeNotValidException {
		super( name, comment, datatype );
		if (datatype != DataType.BOOLEAN && datatype != DataType.BYTE
				&& datatype != DataType.CHAR && datatype != DataType.DOUBLE
				&& datatype != DataType.FLOAT && datatype != DataType.INT
				&& datatype != DataType.LONG && datatype != DataType.SHORT
				&& datatype != DataType.STRING) {
			throw new DatatypeNotValidException(
					"Only primitive datatypes and Strings are allowed as datatypes of a SimpleParameter!");
		}
		this.required = required;
		this.validator = null;
		this.isSet = false;
	}

	/**
	 * Constructor for a {@link SimpleParameter} without
	 * {@link ParameterValidator} but with a default value.
	 * 
	 * @param datatype
	 *            the data type of the parameter value
	 * @param name
	 *            the name of the parameter
	 * @param comment
	 *            a comment on the parameter that tells the user some details
	 *            about it
	 * @param required
	 *            determines if the parameter is required
	 * @param defaultVal
	 *            the default value
	 * 
	 * @throws DatatypeNotValidException if <code>datatype</code> is not in
	 *             the values allowed for a {@link SimpleParameter}
	 * @throws IllegalValueException
	 *             if the default value is not a valid value with respect
	 *             to <code>datatype</code>
	 */
	public SimpleParameter(DataType datatype, String name, String comment,
			boolean required, Object defaultVal) throws DatatypeNotValidException, IllegalValueException {
		this(datatype, name, comment, required);
		setDefault(defaultVal);
	}

	/**
	 * Constructor for a {@link SimpleParameter} with a
	 * {@link ParameterValidator}.
	 * 
	 * @param datatype
	 *            the data type of the parameter value
	 * @param name
	 *            the name of the parameter
	 * @param comment
	 *            a comment on the parameter that tells the user some details
	 *            about it
	 * @param required
	 *            determines if the parameter is required
	 * @param validator
	 *            the validator for the parameter values
	 * 
	 * @throws DatatypeNotValidException
	 *             if <code>datatype</code> is not in the values allowed for a
	 *             {@link SimpleParameter}
	 */
	public SimpleParameter(DataType datatype, String name, String comment,
			boolean required, ParameterValidator validator)
			throws DatatypeNotValidException {
		this(datatype, name, comment, required);
		this.validator = validator;
	}

	/**
	 * Constructor for a {@link SimpleParameter} with validator and default
	 * value.
	 * 
	 * @param datatype
	 *            the data type of the parameter value
	 * @param name
	 *            the name of the parameter
	 * @param comment
	 *            a comment on the parameter that tells the user some details
	 *            about it
	 * @param required
	 *            determines if the parameter is required
	 * @param validator
	 *            the validator for the parameter values
	 * @param defaultVal
	 *            the default value
	 * 
	 * @throws ParameterException
	 *             if either the default value is not a valid value with respect
	 *             to <code>datatype</code> or <code>datatype</code> is not in
	 *             the values allowed for a {@link SimpleParameter}
	 */
	public SimpleParameter(DataType datatype, String name, String comment,
			boolean required, ParameterValidator validator, Object defaultVal)
			throws ParameterException {
		this(datatype, name, comment, required, validator);
		setDefault(defaultVal);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#clone()
	 */
	@Override
	public SimpleParameter clone() throws CloneNotSupportedException {
		SimpleParameter clone = (SimpleParameter) super.clone();
		if (validator != null) {
			clone.validator = validator.clone();
		}
		return clone;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#hasDefaultOrIsSet()
	 */
	@Override
	public boolean hasDefaultOrIsSet() {
		if (isSet) {
			return true;
		} else {
			return getValue() != null;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#isSet()
	 */
	@Override
	public boolean isSet() {
		return isSet;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#isAtomic()
	 */
	@Override
	public boolean isAtomic() {
		return true;
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

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#checkValue(java.lang.Object)
	 */
	@Override
	public boolean checkValue(Object value) {
		if( value == null ) {
			return true;
		}
		boolean res;
		try {
			switch ( datatype ) {
				case BOOLEAN: res = value instanceof Boolean || value instanceof String; break;
				case CHAR: res = value instanceof Character	|| (value instanceof String && ((String) value).length() == 1); break;
				case BYTE: res = value instanceof Byte || (value instanceof String && new Byte((String) value) != null ); break;
				case SHORT: res= value instanceof Short || (value instanceof String && new Short((String) value) != null ); break;
				case INT: res = value instanceof Integer || (value instanceof String && new Integer((String) value) != null ); break;
				case LONG: res = value instanceof Long || (value instanceof String && new Long((String) value) != null ); break;
				case FLOAT: res = value instanceof Float || (value instanceof String && new Float((String) value) != null );break;
				case DOUBLE: res = value instanceof Double || (value instanceof String && new Double((String) value) != null ); break;
				case STRING: res = value instanceof String; break;
				default: res = false;
			}
		} catch ( Exception e ) {
			res = false;
		}
		if( res ) {
			if(validator != null) {
				Object value2 = value;
				if (value instanceof String) {
					//if (((String) value).length() > 0) {
					try {
						switch( datatype ) {
						case BOOLEAN: value2 = new Boolean((String) value); break;
						case BYTE: value2 = new Byte((String) value); break;
						case CHAR: value2 = ((String) value).charAt(0); break;
						case SHORT: value2 = new Short((String) value); break;
						case INT: value2 = new Integer((String) value); break;
						case LONG: value2 = new Long((String) value); break;
						case FLOAT: value2 = new Float((String) value); break;
						case DOUBLE: value2 = new Double((String) value); break;
						//other datatypes are not allowed for this class
						}
					} catch (NumberFormatException e) {
						errorMessage = "Value is not of the expected format.";
						return false;
					}
				}

				res = validator.checkValue(value2);
				errorMessage = validator.getErrorMessage();
			} else {
				errorMessage = "";
			}
		} else {
			errorMessage = "The specified value is no " + datatype + ".";
		}
		return res;
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
	 * @see de.jstacs.parameters.Parameter#setDefault(java.lang.Object)
	 */
	@Override
	public void setDefault(Object defaultValue) throws IllegalValueException {
		if (checkValue(defaultValue)) {
			this.defaultValue = defaultValue;
			setValue(defaultValue);
			isSet = false;
		} else {
			throw new IllegalValueException(name, "Default value ("+defaultValue+") not valid");
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#reset()
	 */
	@Override
	public void reset() {
		if (defaultValue != null) {
			try {
				setDefault(defaultValue);
			} catch (Exception e) {
				value = null;
			}
		} else {
			value = null;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#setValue(java.lang.Object)
	 */
	@Override
	public void setValue(Object value2) throws IllegalValueException {
		if (checkValue(value2)) {
			isSet = true;
			if( value2 == null ) {
				value = null;
				isSet = false;
			} else if (value2 instanceof String && datatype != DataType.STRING) {
				String s = (String) value2;
				// System.out.println("s is: "+s+", datatype: "+datatype);
				try {
					switch (datatype) {
					case BOOLEAN:
						value = new Boolean(s);
						break;
					case CHAR:
						value = s.charAt(0);
						break;
					case BYTE:
						value = new Byte(s);
						break;
					case SHORT:
						value = new Short(s);
						break;
					case INT:
						value = new Integer(s);
						break;
					case LONG:
						value = new Long(s);
						break;
					case FLOAT:
						value = new Float(s);
						break;
					case DOUBLE:
						value = new Double(s);
						break;
					default:
						errorMessage = "Parameter value not of the expected type!";
						throw new IllegalValueException(name,errorMessage);
					}
				} catch (Exception e) {
					this.value = null;
					isSet = false;
					errorMessage = "Value not valid\n" + e.getMessage();
					throw new IllegalValueException(name,"Value not valid\n"
							+ e.getMessage());
				}
			} else {
				this.value = value2;
			}
		} else {
			this.value = null;
			isSet = false;
			throw new IllegalValueException(name,"value not valid: " + value2+"\n"+errorMessage);
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#getValue()
	 */
	@Override
	public Object getValue() {
		return value;
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.AnnotatedEntity#getXMLTag()
	 */
	@Override
	public String getXMLTag() {
		return "simpleParameter";
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#toXML()
	 */
	@Override
	protected void appendFurtherInfos( StringBuffer buf ) {
		super.appendFurtherInfos( buf );
		
		XMLParser.appendObjectWithTags(buf, required, "required");
		XMLParser.appendObjectWithTags(buf, isSet, "isSet");
		XMLParser.appendObjectWithTags(buf, errorMessage, "errorMessage");
		XMLParser.appendObjectWithTags(buf, validator, "validator");
		XMLParser.appendObjectWithTags(buf, defaultValue, "defaultValue");
		XMLParser.appendObjectWithTags(buf, value, "value");
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.parameters.Parameter#extractFurtherInfos(java.lang.StringBuffer)
	 */
	@Override
	protected void extractFurtherInfos( StringBuffer representation ) throws NonParsableException {
		super.extractFurtherInfos( representation );
		
		required = XMLParser.extractObjectForTags(representation, "required", boolean.class );
		errorMessage = XMLParser.parseString( XMLParser.extractObjectForTags(representation, "errorMessage", String.class ) );
		validator = XMLParser.extractObjectForTags( representation, "validator", ParameterValidator.class );
		if( !XMLParser.hasTag(representation, "defaultValue", null, null) ) {
			defaultValue = null;
		} else {
			defaultValue = XMLParser.extractObjectForTags(representation, "defaultValue" );
			if ( defaultValue != null ) {
				try {
					this.setDefault(defaultValue);
				} catch (Exception e) {
					e.printStackTrace();
					throw new NonParsableException(e.getMessage());
				}
			}
		}
		
		String val = XMLParser.extractObjectForTags(representation, "value", String.class );
		if ( val != null ) {
			try {
				this.setValue(val);
			} catch (Exception e) {
				e.printStackTrace();
				throw new NonParsableException(e.getMessage());
			}
		}
		isSet = XMLParser.extractObjectForTags(representation, "isSet", boolean.class );
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#equals(java.lang.Object)
	 */
	@Override
	public boolean equals(Object o2) {
		if (o2 instanceof SimpleParameter) {
			SimpleParameter par2 = (SimpleParameter) o2;
			if (par2.comment.equals(comment) && par2.name.equals(name)
					&& par2.required == required && par2.datatype == datatype) {
				return true;
			} else {
				return false;
			}
		} else {
			return false;
		}
	}

	/**
	 * Class for an {@link Exception} that can be thrown if the provided
	 * <code>int</code>-value that represents a data type is not one of the
	 * values defined in {@link DataType}.
	 * 
	 * @author Jan Grau
	 */
	public static class DatatypeNotValidException extends ParameterException {

		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;

		/**
		 * Creates a new {@link DatatypeNotValidException} with an error
		 * message.
		 * 
		 * @param reason
		 *            the error message
		 */
		public DatatypeNotValidException(String reason) {
			super("The datatype is not valid for this type of parameter: "
					+ reason);
		}

	}

	/**
	 * This exception is thrown if a parameter is not valid.
	 * 
	 * @author Jan Grau
	 */
	public static class IllegalValueException extends ParameterException {

		/**
		 * 
		 */
		private static final long serialVersionUID = 1L;

		/**
		 * Creates a new {@link IllegalValueException} with the reason of the
		 * exception <code>reason</code> as error message.
		 * 
		 * @param reason
		 *            the reason to give as error message
		 */
		public IllegalValueException(String name, String reason) {
			super(name, "Parameter not permitted: " + reason);
		}

	}

	/**
	 * Returns the {@link ParameterValidator} used in this
	 * {@link SimpleParameter}. This may be <code>null</code>.
	 * 
	 * @return the validator used in this {@link SimpleParameter}
	 */
	public ParameterValidator getValidator() {
		return validator;
	}

	/**
	 * Sets a new {@link ParameterValidator} for this {@link SimpleParameter}.
	 * 
	 * @param validator
	 *            the new parameter validator
	 */
	public void setValidator(ParameterValidator validator) {
		this.validator = validator;
	}
	
	/*
	 * (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	public String toString() {
		return name + " (" + comment
			+ (validator!=null && validator instanceof NumberValidator?", " + validator.toString():"")
			+ (defaultValue!=null?", default = " + defaultValue.toString():"")
			+ (required ? "" : ", OPTIONAL" )
			+ ")\t= " + value;
	}

	@Override
	public void toGalaxy( String namePrefix, String configPrefix, int depth, StringBuffer descBuffer, StringBuffer configBuffer, boolean addLine, int indentation ) throws Exception {
		namePrefix = namePrefix+"_"+GalaxyAdaptor.getLegalName( getName() );
		int nextIndentation = XMLParser.nextIndentation(indentation);
		StringBuffer buf = new StringBuffer();
		if(validator != null && validator instanceof GalaxyConvertible){
			((GalaxyConvertible)validator).toGalaxy( namePrefix+"_valid", null, depth, buf, null, false, nextIndentation );
		}
		
		String line = "";
		if(addLine){
			line = "&lt;hr /&gt;";
		}
		
		XMLParser.addTagsAndAttributes( buf, "param", "type=\""+dataTypeToGalaxy()+"\""+(datatype == DataType.STRING ? " size=\"40\"" : "")+" name=\""+namePrefix+"\" label=\""+line+getName()+"\" help=\""+XMLParser.escape(getComment())+"\" "+(datatype == DataType.BOOLEAN ? "checked" : "value")+"=\""+(defaultValue == null ? "" : (datatype == DataType.BOOLEAN ? (defaultValue.equals( true ) ? "True" : "False") : defaultValue) )+"\" optional=\""+(!isRequired())+"\"", indentation );
		descBuffer.append( buf );
		
		buf = new StringBuffer();
		buf.append( "${"+configPrefix+namePrefix+"}" );
		XMLParser.addTags( buf, namePrefix );
		configBuffer.append( buf );		
	}
	
	

	/**
	 * Returns the Galaxy identifier for the {@link DataType} of this parameter
	 * @return the Galaxy identifier of the data type
	 */
	protected String dataTypeToGalaxy() {
		
		switch(datatype){
			case LONG:
			case INT:
			case SHORT:
			case BYTE:
				return "integer";
			case FLOAT:
			case DOUBLE:
				return "float";
			case BOOLEAN:
				return "boolean";
			case CHAR:
			case STRING:
				return "text";
			default:
				return "text";		
		}
		
	}

	@Override
	public void fromGalaxy( String namePrefix, StringBuffer command ) throws Exception {
		namePrefix = namePrefix+"_"+GalaxyAdaptor.getLegalName( getName() );
		StringBuffer value = XMLParser.extractForTag( command, namePrefix );
		
		if( value != null ) {
			String	val = unescape( value.toString() );
			if( isRequired() || val.length()>0 || checkValue(val) ) this.setValue(val); 
		} else {
			if( isRequired() && !hasDefaultOrIsSet() ) {
				throw new NullPointerException( getName()+" "+command+" "+namePrefix );
			}
		}
	}
	
	private static String[][] table = {
	                                   { "&", "__amp__" },
	                                   { "\"", "__quot__" },
	                                   { "'", "__sq__"}, //"__apos__" },
	                                   { ">", "__gt__" },
	                                   { "<", "__lt__" },
	                                   { "\n", "(__cn__|__cr__)+"},
	                                   { "[", "__ob__" },
	                                   { "]", "__cb__" }
	};

	private static String unescape( String original ) {
		//System.out.println("before: "+original);
		if( original != null ) {
			for( int i = 0; i < table.length; i++ ) {
				original = original.replaceAll( table[i][1], table[i][0] );
			}
		}
		//System.out.println("after: "+original);
		return original;
	}

	@Override
	public void toGalaxyTest(String namePrefix, int depth, StringBuffer testBuffer, int indentation) throws Exception {
		if( isSet ) {
			XMLParser.insertIndentation(testBuffer, indentation, testBuffer.length());
			testBuffer.append("<param name=\"" + namePrefix+"_"+GalaxyAdaptor.getLegalName( getName() ) + "\" value=\"" + getValue().toString() + "\" />\n");
		}
	}
}
