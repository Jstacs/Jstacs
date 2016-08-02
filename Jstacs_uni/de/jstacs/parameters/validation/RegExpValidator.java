package de.jstacs.parameters.validation;

import java.util.regex.Pattern;

import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.GalaxyConvertible;

public class RegExpValidator implements ParameterValidator, GalaxyConvertible {

	private String regExp;
	private String errorMessage;
	
	public RegExpValidator(String regExp){
		this.regExp = regExp;
	}
	
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
			StringBuffer configBuffer, boolean addLine) throws Exception {
		descBuffer.append( "<validator type=\"regex\" message=\"Value does not match "+regExp+".\">^"+regExp+"$</validator>" );
		
	}

	@Override
	public void fromGalaxy(String namePrefix, StringBuffer command) throws Exception {
		
	}

}
