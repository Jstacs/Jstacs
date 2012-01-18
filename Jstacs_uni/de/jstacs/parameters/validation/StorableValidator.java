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

package de.jstacs.parameters.validation;

import de.jstacs.Storable;
import de.jstacs.classifiers.AbstractClassifier;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.FileParameter;
import de.jstacs.sequenceScores.statisticalModels.trainable.AbstractTrainSM;
import de.jstacs.sequenceScores.statisticalModels.trainable.TrainableStatisticalModel;

/**
 * Class for a validator that validates instances and XML representations for
 * the correct class types (e.g. {@link de.jstacs.sequenceScores.statisticalModels.trainable.AbstractTrainSM}).
 * 
 * @author Jan Grau
 */
public class StorableValidator implements ParameterValidator {

	/**
	 * The error message, <code>null</code> if no error occurred
	 */
	private String errorMessage;
	/**
	 * The class of the object
	 */
	private Class clazz;
	/**
	 * <code>true</code> if the model or classifier must be trained
	 */
	private boolean trained;

	/**
	 * Creates a new {@link StorableValidator} for a subclass of
	 * {@link AbstractTrainSM} or {@link AbstractClassifier}. This constructor may
	 * not be used on other subclasses of {@link Storable}.
	 * 
	 * @param clazz
	 *            the class
	 * @param trained
	 *            <code>true</code> if the model or classifier must be trained
	 * @throws Exception
	 *             if <code>clazz</code> is not of the expected type
	 */
	public StorableValidator(Class<? extends Storable> clazz, boolean trained)
			throws Exception {
		if (AbstractTrainSM.class.isAssignableFrom(clazz)
				|| AbstractClassifier.class.isAssignableFrom(clazz)) {
			this.clazz = clazz;
		} else {
			throw new Exception(
					"Class is not a subtype of AbstractClassifier or AbstractTrainSM.");
		}
		this.trained = trained;
	}

	/**
	 * Creates a new {@link StorableValidator} for a subclass of
	 * {@link Storable}.
	 * 
	 * @param clazz
	 *            the class
	 * 
	 * @throws Exception
	 *             if <code>clazz</code> is not of the expected type
	 */
	public StorableValidator(Class<? extends Storable> clazz) throws Exception {
		if (Storable.class.isAssignableFrom(clazz)) {
			this.clazz = clazz;
		} else {
			throw new Exception("Class is not a subtype of Storable");
		}
		this.trained = false;
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Constructs an {@link StorableValidator} from its XML representation.
	 * 
	 * @param buf
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if <code>buf</code> could not be parsed
	 */
	public StorableValidator(StringBuffer buf) throws NonParsableException {
		fromXML(buf);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#clone()
	 */
	@Override
	public StorableValidator clone() throws CloneNotSupportedException {
		try {
			StorableValidator clone = null;
			if (AbstractTrainSM.class.isAssignableFrom(clazz)
					|| AbstractClassifier.class.isAssignableFrom(clazz)) {
				clone = new StorableValidator(
						(Class<? extends Storable>) clazz, trained);
			} else {
				clone = new StorableValidator((Class<? extends Storable>) clazz);
			}
			clone.errorMessage = errorMessage;
			return clone;
		} catch (Exception e) {
			throw new CloneNotSupportedException(e.getCause().getMessage());
		}
	}

	/**
	 * Checks the value of <code>value</code>. Allowed types of
	 * <code>value</code> are {@link AbstractTrainSM}, {@link AbstractClassifier},
	 * {@link de.jstacs.parameters.FileParameter.FileRepresentation}, {@link String}, and {@link StringBuffer}. In
	 * all cases where an XML representation is given as <code>value</code>, it
	 * must be surrounded by &lt;object&gt;-tags and these tags must contain a
	 * &lt;className&gt;-element that contains the name of the class of the
	 * represented instance.
	 * 
	 * @param value
	 *            the {@link Object} to be checked
	 * 
	 * @return <code>true</code> if <code>value</code> is valid and
	 *         <code>false</code> otherwise
	 */
	@SuppressWarnings("unchecked")
	public boolean checkValue(Object value) {
		if (AbstractTrainSM.class.isAssignableFrom(clazz)
				&& value instanceof AbstractTrainSM) {
			if (!trained || ((TrainableStatisticalModel) value).isInitialized()) {
				errorMessage = null;
				return true;
			} else {
				errorMessage = "The model must be trained.";
				return false;
			}
		} else if (AbstractClassifier.class.isAssignableFrom(clazz)
				&& value instanceof AbstractClassifier) {
			if (!trained || ((AbstractClassifier) value).isInitialized()) {
				errorMessage = null;
				return true;
			} else {
				errorMessage = "The classifier must be trained.";
				return false;
			}
		} else {

			StringBuffer buf = null;
			if (value instanceof FileParameter.FileRepresentation) {
				String content = ((FileParameter.FileRepresentation) value)
						.getContent();
				buf = new StringBuffer(content);
			} else if (value instanceof String) {
				buf = new StringBuffer((String) value);
			} else if (value instanceof StringBuffer) {
				buf = (StringBuffer) value;
			}

			try {
				buf = XMLParser.extractForTag(buf, "object");
				String className = XMLParser.extractObjectForTags(buf, "className", String.class );
				Class c = Class.forName(className);
				if (clazz.isAssignableFrom(c)) {
					boolean modelTrained = false;
					if (AbstractTrainSM.class.isAssignableFrom(c)) {
						TrainableStatisticalModel model = (TrainableStatisticalModel) c.getConstructor(
								new Class[] { StringBuffer.class })
								.newInstance(buf);
						modelTrained = model.isInitialized();
					} else if (AbstractClassifier.class.isAssignableFrom(c)) {
						AbstractClassifier classifier = (AbstractClassifier) c
								.getConstructor(
										new Class[] { StringBuffer.class })
								.newInstance(buf);
						modelTrained = classifier.isInitialized();
					} else {
						c.getConstructor(new Class[] { StringBuffer.class })
								.newInstance(buf);
					}
					if (!trained || modelTrained) {
						errorMessage = null;
						return true;
					} else {
						errorMessage = "The classifier must be trained.";
						return false;
					}
				} else {
					errorMessage = "File content was not of the correct class "
							+ clazz.getName() + " but of " + className;
					return false;
				}
			} catch (Exception e) {
				e.printStackTrace();
				errorMessage = "The file could not be parsed: " + e.getCause();
				return false;
			}
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.validation.ParameterValidator#getErrorMessage()
	 */
	public String getErrorMessage() {
		return errorMessage;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.Storable#toXML()
	 */
	public StringBuffer toXML() {
		StringBuffer buf = new StringBuffer();
		XMLParser.appendObjectWithTags(buf, errorMessage, "errorMessage");
		XMLParser.appendObjectWithTags(buf, trained, "trained");
		XMLParser.appendObjectWithTags(buf, clazz.getName(), "class");
		XMLParser.addTags(buf, "TrainedValidator");

		return buf;
	}

	/**
	 * Parses a {@link StorableValidator} from the XML representation as
	 * returned by {@link StorableValidator#toXML()}.
	 * 
	 * @param representation
	 *            the XML representation
	 * 
	 * @throws NonParsableException
	 *             if the XML code could not be parsed
	 */
	public void fromXML(StringBuffer representation)
			throws NonParsableException {
		representation = XMLParser.extractForTag(representation,
				"TrainedValidator");
		errorMessage = XMLParser.extractObjectForTags(representation, "errorMessage", String.class );
		trained = XMLParser.extractObjectForTags(representation, "trained", boolean.class );
		try {
			clazz = Class.forName(XMLParser.extractObjectForTags(representation, "class", String.class ));
		} catch (Exception e) {
			throw new NonParsableException(e.getMessage());
		}
	}

}
