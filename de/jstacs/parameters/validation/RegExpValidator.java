package de.jstacs.parameters.validation;

import java.util.regex.Pattern;

import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.GalaxyConvertible;

/**
 * {@link ParameterValidator} that checks if a given input {@link String} matches a regular expression.
 * 
 * @author Jan Grau, Jens Keilwagen
 *
 */
public class RegExpValidator implements ParameterValidator, GalaxyConvertible {

	private String regExp;
	private String errorMessage;
	
	/**
	 * Creates a validator for the given regular expression.
	 * @param regExp the regular expression
	 * @see Pattern
	 */
	public RegExpValidator(String regExp){
		this.regExp = regExp;
	}
	
	/**
	 * Creates a {@link RegExpValidator} from its XML representation.
	 * @param xml the XML representation
	 * @throws NonParsableException if the XML representation could not be parsed
	 */
	public RegExpValidator(StringBuffer xml) throws NonParsableException{
		xml = XMLParser.extractForTag(xml, "RegExpValidator");
		regExp = (String) XMLParser.extractObjectForTags(xml, "regExp");
	}
	
	@Override
	public StringBuffer toXML() {
		StringBuffer xml = new StringBuffer();
		XMLParser.appendObjectWithTags(xml, regExp, "regExp");
		XMLParser.addTags(xml, "RegExpValidator");
		return xml;
	}

	@Override
	public boolean checkValue(Object value) {
		if(value == null){
			errorMessage = "Value is null.";
			return false;
		}else if(!(value instanceof String)){
			errorMessage = "Value is no String.";
			return false;
		}else{
			String str = (String) value;
			if( Pattern.matches(regExp, str) ){
				errorMessage = null;
				return true;
			}else{
				errorMessage = "String does not match "+regExp;
				return false;
			}
		}
	}

	@Override
	public String getErrorMessage() {
		return errorMessage;
	}

	@Override
	public RegExpValidator clone() throws CloneNotSupportedException {
		RegExpValidator val = (RegExpValidator) super.clone();
		return val;
	}

	@Override
	public void toGalaxy(String namePrefix, String configPrefix, int depth, StringBuffer descBuffer,
			StringBuffer configBuffer, boolean addLine, int indentation) throws Exception {
		XMLParser.addIndentation(descBuffer, indentation);
		descBuffer.append( "<validator type=\"regex\" message=\"Value does not match "+XMLParser.escape(regExp)+".\">^"+XMLParser.escape(regExp)+"$</validator>\n" );
	}

	@Override
	public void fromGalaxy(String namePrefix, StringBuffer command) throws Exception {
		
	}

}
