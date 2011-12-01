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
import de.jstacs.Storable;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;
import de.jstacs.parameters.validation.ParameterValidator;
import de.jstacs.utils.galaxy.GalaxyAdaptor;

/**
 * Class for a {@link Parameter} that represents a local file.
 * 
 * @author Jan Grau
 * 
 */
public class FileParameter extends Parameter implements GalaxyConvertible {

	/**
	 * The name of the parameter
	 */
	private String name;
	/**
	 * The comment on the parameter
	 */
	private String comment;
	/**
	 * <code>true</code> if the parameter is required
	 */
	private boolean required;
	/**
	 * The MIME-type of accepted files
	 */
	private String mime;
	/**
	 * The file
	 */
	private FileRepresentation value;
	/**
	 * The default value
	 */
	private FileRepresentation defaultValue;
	/**
	 * <code>true</code> if a file is set as value
	 */
	private boolean isSet;
	/**
	 * The error message, <code>null</code> if no error occurred
	 */
	private String errorMessage;

	/**
	 * The parameter validator
	 */
	private ParameterValidator valid;

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#clone()
	 */
	@Override
	public FileParameter clone() throws CloneNotSupportedException {
		FileParameter clone = (FileParameter) super.clone();
		clone.value = value == null ? null : value.clone();
		clone.defaultValue = defaultValue == null ? null : defaultValue.clone();
		return clone;
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Restores a {@link FileParameter} from an XML representation.
	 * 
	 * @param buf
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the XML code could not be parsed
	 * 
	 * @see FileParameter#fromXML(StringBuffer)
	 */
	public FileParameter(StringBuffer buf) throws NonParsableException {
		fromXML(buf);
	}

	/**
	 * Creates a {@link FileParameter}.
	 * 
	 * @param name
	 *            the name of the parameter
	 * @param comment
	 *            a comment on the parameter
	 * @param mime
	 *            the MIME-type of allowed files
	 * @param required
	 *            <code>true</code> if this {@link FileParameter} is required to
	 *            continue, <code>false</code> otherwise
	 */
	public FileParameter(String name, String comment, String mime,
			boolean required) {
		this.name = name;
		this.comment = comment;
		this.mime = mime;
		this.required = required;
	}

	/**
	 * Constructs a {@link FileParameter}.
	 * 
	 * @param name
	 *            the name of the parameter
	 * @param comment
	 *            a comment on the parameter
	 * @param mime
	 *            the MIME-type of allowed files
	 * @param required
	 *            <code>true</code> if this {@link FileParameter} is required
	 * @param validator
	 *            a validator that validates e.g. the contents of the file
	 */
	public FileParameter(String name, String comment, String mime,
			boolean required, ParameterValidator validator) {
		this(name, comment, mime, required);
		this.valid = validator;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#getName()
	 */
	@Override
	public String getName() {
		return name;
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
	 * @see de.jstacs.parameters.Parameter#getDatatype()
	 */
	@Override
	public DataType getDatatype() {
		return DataType.FILE;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#getComment()
	 */
	@Override
	public String getComment() {
		return comment;
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
	 * Resets the {@link FileParameter} to its original state.
	 */
	@Override
	public void reset() {
		this.value = null;
		this.isSet = false;
		this.errorMessage = null;
	}

	/**
	 * Returns the content of the file.
	 * 
	 * @return the content of the file
	 */
	public FileRepresentation getFileContents() {
		return value;
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
	 * @see de.jstacs.parameters.Parameter#checkValue(java.lang.Object)
	 */
	@Override
	public boolean checkValue(Object value) {
		if (valid != null) {
			if (valid.checkValue(value)) {
				errorMessage = null;
				return true;
			} else {
				errorMessage = valid.getErrorMessage();
				return false;
			}
		} else if (value != null && value instanceof FileRepresentation) {
			FileRepresentation f = (FileRepresentation) value;
			if (f.getFilename() != null && f.getFilename().length() != 0
					&& f.getContent() != null && f.getContent().length() != 0) {
				errorMessage = null;
				return true;
			} else {
				errorMessage = "No file specified or file is empty.";
				return false;
			}
		} else {
			errorMessage = "Value is no file or null.";
			return false;
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#setDefault(java.lang.Object)
	 */
	@Override
	public void setDefault(Object defaultValue) throws IllegalValueException {
		if (checkValue(defaultValue)) {
			this.defaultValue = (FileRepresentation) defaultValue;
			setValue(defaultValue);
			isSet = false;
		} else {
			throw new IllegalValueException(errorMessage);
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#simplify()
	 */
	@Override
	public void simplify() {

	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#setValue(java.lang.Object)
	 */
	@Override
	public void setValue(Object value) throws IllegalValueException {
		if (!checkValue(value)) {
			throw new IllegalValueException(errorMessage);
		}
		this.value = (FileRepresentation) value;
		this.isSet = true;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#getValue()
	 */
	@Override
	public Object getValue() {
		if (value == null) {
			return null;
		} else {
			return value.getFilename();
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#hasDefaultOrIsSet()
	 */
	@Override
	public boolean hasDefaultOrIsSet() {
		return isSet();
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
	 * @see de.jstacs.parameters.Parameter#toXML()
	 */
	@Override
	public StringBuffer toXML() {
		StringBuffer buf = super.toXML();
		XMLParser.addTags(buf, "superParameter");
		XMLParser.appendObjectWithTags(buf, name, "name");
		XMLParser.appendObjectWithTags(buf, comment, "comment");
		XMLParser.appendObjectWithTags(buf, mime, "mime");
		XMLParser.appendObjectWithTags(buf, required, "required");
		XMLParser.appendObjectWithTags(buf, isSet, "isSet");
		XMLParser.appendObjectWithTags(buf, errorMessage, "errorMessage");
		
		XMLParser.appendObjectWithTags(buf, value,"value");
		XMLParser.appendObjectWithTags(buf, valid, "validator");
		/*
		if (value == null) {
			XMLParser.appendObjectWithTags(buf, "null", "value");
		} else {
			XMLParser.appendObjectWithTags(buf, value.toXML().toString(),
					"value");
		}
		if (valid == null) {
			XMLParser.appendObjectWithTags(buf, "null", "validator");
		} else {
			StringBuffer buf2 = new StringBuffer();
			XMLParser.appendObjectWithTags(buf2, valid.getClass().getName(),
					"className");
			buf2.append(valid.toXML());
			XMLParser.appendObjectWithTags(buf, buf2.toString(), "validator");
		}
		*/
		XMLParser.addTags(buf, "fileParameter");

		return buf;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#fromXML(java.lang.StringBuffer)
	 */
	@Override
	protected void fromXML(StringBuffer representation)
			throws NonParsableException {
		StringBuffer buf = XMLParser.extractForTag(representation,"fileParameter");
		super.fromXML(XMLParser.extractForTag(buf,"superParameter"));
		name = XMLParser.extractObjectForTags(buf, "name", String.class );// TODO XMLP14CONV This and (possibly) the following lines have been converted automatically
		comment = XMLParser.extractObjectForTags(buf, "comment", String.class );
		mime = XMLParser.extractObjectForTags(buf, "mime", String.class );
		required = XMLParser.extractObjectForTags(buf, "required", boolean.class );
		isSet = XMLParser.extractObjectForTags(buf, "isSet", boolean.class );
		errorMessage = XMLParser.parseString( XMLParser.extractObjectForTags(buf, "errorMessage", String.class ) );
		
		value = XMLParser.extractObjectForTags(buf, "value", FileRepresentation.class );
		valid = XMLParser.extractObjectForTags(buf, "validator", ParameterValidator.class );
		/*
		String val = XMLParser.extractObjectForTags(buf, "value", String.class );
		if (val == null ) {//XXX
			value = null;
		} else {
			value = new FileRepresentation(new StringBuffer(val));
		}
		System.out.println("hier");
		System.exit( 1 );
		val = XMLParser.extractObjectForTags(buf, "validator", String.class );
		if (val.equals("null")) {
			valid = null;
		} else {
			StringBuffer buf2 = new StringBuffer(val);
			String className = XMLParser.extractObjectForTags(buf2, "className", String.class );// TODO XMLP14CONV This and (possibly) the following lines have been converted automatically
			try {
				valid = (ParameterValidator) Class.forName(className)
						.getConstructor(new Class[] { StringBuffer.class })
						.newInstance(buf2);
			} catch (Exception e) {
				throw new NonParsableException(e.getMessage());
			}
		}
		*/
	}

	/**
	 * Returns the MIME-type of the allowed files.
	 * 
	 * @return the MIME-type of the allowed files
	 */
	public String getAcceptedMimeType() {
		return mime;
	}
	
	@Override
	public void toGalaxy( String namePrefix, String configPrefix, int depth, StringBuffer descBuffer, StringBuffer configBuffer  ) {
		namePrefix = namePrefix+"_"+GalaxyAdaptor.getLegalName( getName() );
		StringBuffer buf = new StringBuffer();
		XMLParser.addTagsAndAttributes( buf, "param", "type=\"data\" name=\""+namePrefix+"\" label=\""+getName()+"\" help=\""+getComment()+"\" value=\""+(defaultValue == null ? "" : defaultValue)+"\" optional=\""+(!isRequired())+"\"" );
		descBuffer.append( buf );
		buf = new StringBuffer();
		buf.append( "${"+configPrefix+namePrefix+"}" );
		XMLParser.addTags( buf, namePrefix );
		configBuffer.append( buf );
	}

	@Override
	public void fromGalaxy( String namePrefix, StringBuffer command ) throws Exception {
		namePrefix = namePrefix+"_"+GalaxyAdaptor.getLegalName( getName() );
		String val = XMLParser.extractForTag( command, namePrefix ).toString();
		this.value = new FileRepresentation( val );
		this.isSet = true;
	}
	
	/**
	 * Class that represents a file.
	 * 
	 * @author Jan Grau
	 */
	public static class FileRepresentation implements Storable, Cloneable {

		/**
		 * The name of the file
		 */
		private String filename;
		/**
		 * The contents of the file
		 */
		private String content;

		/**
		 * Creates a {@link FileRepresentation} out of the filename and the
		 * file's contents.
		 * 
		 * @param filename
		 *            the name of the file
		 * @param content
		 *            the contents of the file
		 */
		public FileRepresentation(String filename, String content) {
			this.filename = filename;
			this.content = content;
		}
		
		/**
		 * Creates a new {@link FileRepresentation} from a filename.
		 * @param filename the filename
		 */
		public FileRepresentation(String filename){
			this.filename = filename;
		}

		/**
		 * The standard constructor for the interface {@link de.jstacs.Storable}
		 * . Restores the {@link FileRepresentation} from an XML representation.
		 * 
		 * @param buf
		 *            the XML representation as {@link StringBuffer}
		 * 
		 * @throws NonParsableException
		 *             if the {@link StringBuffer} could not be parsed
		 * 
		 * @see #fromXML(StringBuffer)
		 */
		public FileRepresentation(StringBuffer buf) throws NonParsableException {
			fromXML(buf);
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see java.lang.Object#clone()
		 */
		@Override
		public FileRepresentation clone() throws CloneNotSupportedException {
			return (FileRepresentation) super.clone();
		}

		/**
		 * Returns the filename.
		 * 
		 * @return the name of the file
		 */
		public String getFilename() {
			return filename;
		}

		/**
		 * Returns the content of the file.
		 * 
		 * @return the content of the file
		 */
		public String getContent() {
			return content;
		}

		/*
		 * (non-Javadoc)
		 * 
		 * @see de.jstacs.Storable#toXML()
		 */
		public StringBuffer toXML() {
			StringBuffer buf = new StringBuffer();
			XMLParser.appendObjectWithTags(buf, filename, "filename");
			XMLParser.appendObjectWithTags(buf, content, "content");
			XMLParser.addTags(buf, "fileRepresentation");

			return buf;
		}

		private void fromXML(StringBuffer representation)
				throws NonParsableException {
			representation = XMLParser.extractForTag(representation,
					"fileRepresentation");
			filename = XMLParser
					.extractObjectForTags(representation, "filename", String.class );// TODO XMLP14CONV This and (possibly) the following lines have been converted automatically
			content = XMLParser.extractObjectForTags(representation, "content", String.class );
		}

	}

}
