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
import de.jstacs.parameters.SimpleParameter.DatatypeNotValidException;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;
import de.jstacs.parameters.validation.NumberValidator;

/**
 * Class for a {@link ParameterSet} that consists of a length-{@link Parameter}
 * that defines the length of the array and an array of
 * {@link ParameterSetContainer}s of this length. Each
 * {@link ParameterSetContainer} in the array holds a {@link ParameterSet} as
 * defined by a template. The template is cloned in order to obtain the
 * specified number of elements in the array. The class takes care of the
 * consistency of the length-{@link Parameter} and the length of the array.
 * 
 * @author Jan Grau
 * 
 */
public class ArrayParameterSet extends ExpandableParameterSet {

	private int numberChanged;
	private String lengthName;
	private String lengthComment;
	private NumberValidator<Integer> allowedLengths;

	/**
	 * Creates a new {@link ArrayParameterSet} from a {@link Class} that can be
	 * instantiated using this {@link ArrayParameterSet} and templates for the
	 * {@link ParameterSet} in each element of the array, the name and the
	 * comment that are displayed for the {@link ParameterSetContainer}s
	 * enclosing the {@link ParameterSet}s.
	 * 
	 * @param template
	 *            the template of the {@link ParameterSet}
	 * @param nameTemplate
	 *            the name-template
	 * @param commentTemplate
	 *            the comment-template
	 * @param lengthName
	 *            the name of the length-{@link Parameter}
	 * @param lengthComment
	 *            the comment of the length-{@link Parameter}
	 * @param allowedLengths
	 *            a {@link NumberValidator} to set a lower and upper bound on
	 *            the number of elements in the array
	 * @throws CloneNotSupportedException if the templace could not be cloned
	 * @throws IllegalValueException if the values of <code>allowedLengths</code> are not allowed values 
	 * 		of the length parameter
	 * 
	 */
	public ArrayParameterSet(ParameterSet template, String nameTemplate,
			String commentTemplate, String lengthName, String lengthComment,
			NumberValidator<Integer> allowedLengths) throws CloneNotSupportedException, IllegalValueException {
		super(template, nameTemplate, commentTemplate);
		this.lengthName = lengthName;
		this.lengthComment = lengthComment;
		this.allowedLengths = allowedLengths;

		SimpleParameter length;
		try {
			length = new SimpleParameter(DataType.INT,
					lengthName, lengthComment, true, allowedLengths);
		} catch ( DatatypeNotValidException doesnothappen ) {
			throw new RuntimeException( doesnothappen );
		}
		length.setDefault(allowedLengths.getLowerBound());
		length.setValue(allowedLengths.getLowerBound());
		length.setRangeable(false);
		this.parameters.add(0,length);
		//this.addParameterToSet();//TODO
		for(int i=1;i<(Integer)length.getValue();i++){
			this.addParameterToSet();
		}
		numberChanged = 1;

	}

	/**
	 * Creates a new {@link ArrayParameterSet} from a {@link Class} that can be
	 * instantiated using this {@link ArrayParameterSet} and templates for the
	 * {@link ParameterSet} in each element of the array, the name and the
	 * comment that are displayed for the {@link ParameterSetContainer}s
	 * enclosing the {@link ParameterSet}s.
	 * 
	 * @param template
	 *            the template of the {@link ParameterSet}
	 * @param nameTemplate
	 *            the name-template
	 * @param commentTemplate
	 *            the comment-template
	 * @throws CloneNotSupportedException if the template could not be cloned
	 * @throws IllegalValueException if the values of <code>allowedLengths</code> are not allowed values 
	 * 		of the length parameter
	 */
	public ArrayParameterSet(ParameterSet template, String nameTemplate,
			String commentTemplate) throws IllegalValueException, CloneNotSupportedException {
		this(template, nameTemplate, commentTemplate, "Length",
				"The length of the array.", new NumberValidator<Integer>(1,
						Integer.MAX_VALUE));
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link ArrayParameterSet} from its XML representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link StringBuffer} <code>representation</code> could
	 *             not be parsed
	 */
	public ArrayParameterSet(StringBuffer representation)
			throws NonParsableException {
		super(representation);
	}
	
	//TODO length may be set directly in parameter. Not handled so far.
	public void setLength(int length) throws IllegalValueException, CloneNotSupportedException{
		this.parameters.get(0).setValue((Integer)length);
		int numBefore = this.parameters.size()-1;
		if(numBefore < length){
			for(int i=numBefore;i<length;i++){
				this.addParameterToSet();
			}
		}else if(numBefore > length){
			for(int i=numBefore;i>length;i--){
				this.removeParameterFromSet();
			}
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.ParameterSet#getNumberOfParameters()
	 */
	@Override
	public int getNumberOfParameters() {
		return super.getNumberOfParameters();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.ParameterSet#getParameterAt(int)
	 */
	@Override
	public Parameter getParameterAt(int i) {
		return super.getParameterAt(i);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.ParameterSet#hasDefaultOrIsSet()
	 */
	@Override
	public boolean hasDefaultOrIsSet() {
		if (numberChanged != 0) {
			return false;
		} else {
			return super.hasDefaultOrIsSet();
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.ExpandableParameterSet#fromXML(java.lang.StringBuffer)
	 */
	@Override
	protected void fromXML(StringBuffer representation)
			throws NonParsableException {
		representation = XMLParser.extractForTag(representation,
				"arrayParameterSet");

		super.fromXML(XMLParser.extractForTag(representation,
				"superRepresentation"));
		lengthName = XMLParser
				.extractObjectForTags(representation, "lengthName", String.class );
		lengthComment = XMLParser.extractObjectForTags(representation, "lengthComment", String.class );
		allowedLengths = XMLParser.extractObjectForTags(representation, "allowedLengths", NumberValidator.class );
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.ExpandableParameterSet#toXML()
	 */
	@Override
	public StringBuffer toXML() {
		StringBuffer supBuf = super.toXML();
		XMLParser.addTags(supBuf, "superRepresentation");
		XMLParser.appendObjectWithTags(supBuf, lengthName, "lengthName");
		XMLParser.appendObjectWithTags(supBuf, lengthComment, "lengthComment");
		XMLParser.appendObjectWithTags(supBuf, allowedLengths,
				"allowedLengths");
		XMLParser.addTags(supBuf, "arrayParameterSet");
		return supBuf;

	}

}
