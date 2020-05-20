package de.jstacs.parameters.validation;

import java.io.File;

import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.FileParameter.FileRepresentation;

/**
 * A simple validator allowing to check whether a file used in a FileParameter exists.
 * 
 * @author Jens Keilwagen
 * 
 * @see FileParameter
 */
public class FileExistsValidator implements ParameterValidator {

	private String errorMsg;
	
	/**
	 * Standard constructor.
	 */
	public FileExistsValidator() {}
	
	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link FileExistsValidator} from its XML representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if <code>representation</code> could not be parsed
	 */
	public FileExistsValidator(StringBuffer xml) throws NonParsableException {
		StringBuffer part = XMLParser.extractForTag(xml, "FileExistsValidator");
		errorMsg=XMLParser.extractObjectForTags(part, "errorMsg", String.class);
	}
	
	@Override
	public StringBuffer toXML() {
		StringBuffer xml = new StringBuffer();
		XMLParser.appendObjectWithTags(xml, errorMsg, "errorMsg");
		XMLParser.addTags(xml, "FileExistsValidator");
		return xml;
	}

	@Override
	public boolean checkValue(Object value) {
		if (value == null) {
			errorMsg = "Value is null.";
			return false;
		}
		if( value instanceof FileRepresentation ) {
			FileRepresentation fr = (FileRepresentation) value;
			boolean ex = new File(fr.getFilename()).exists();
			if( !ex ) {
				errorMsg="File " + fr.getFilename() + " does not exist";
			}
			return ex;
		} else {
			errorMsg="The object is no FileRepresentation.";
			return false;
		}
	}

	@Override
	public String getErrorMessage() {
		return errorMsg;
	}

	@Override
	public FileExistsValidator clone() throws CloneNotSupportedException {
		return (FileExistsValidator) super.clone();
	}
}