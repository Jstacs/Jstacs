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

import java.util.StringTokenizer;

import de.jstacs.DataType;
import de.jstacs.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;
import de.jstacs.utils.galaxy.GalaxyAdaptor;

/**
 * Class for a parameter wrapper that allows {@link SimpleParameter}s to be set
 * to a set of values.<br>
 * These values may be given either as a list of values separated by spaces, as
 * a range between a first and a last value with a given number of steps between
 * these values, or a single value. In the latter case a {@link RangeParameter}
 * works just like a plain {@link SimpleParameter}. {@link SimpleParameter}s
 * that are used to construct a {@link RangeParameter} must be rangeable, i.e.
 * their method {@link SimpleParameter#isRangeable()} must return
 * <code>true</code>, or they will be rejected.
 * 
 * @author Jan Grau, Jens Keilwagen
 */
public class RangeParameter extends Parameter implements RangeIterator, GalaxyConvertible {

	/**
	 * This <code>enum</code> defines possible scales, if {@link RangeType} is
	 * {@link RangeType#RANGE}.
	 * 
	 * @author Jan Grau
	 */
	public enum Scale {
		/**
		 * The range is sampled in a logarithmic scale.
		 */
		LOGSCALE,
		/**
		 * The range is sampled in an inverse (finer around the last value)
		 * logarithmic scale.
		 */
		INVERSELOGSCALE,
		/**
		 * The range is sampled in a linear scale.
		 */
		LINSCALE,
		/**
		 * The range is sampled in an exponential scale.
		 */
		EXPSCALE,
		/**
		 * The range is sampled in a logarithmic scale (radix=10).
		 */
		LOGSCALE2,
		/**
		 * The range is sampled in an inverse (finer around the last value)
		 * logarithmic scale (radix=2).
		 */
		INVERSELOGSCALE2,
		/**
		 * The range is sampled in a logarithmic scale (radix=10).
		 */
		LOGSCALE10,
		/**
		 * The range is sampled in an inverse (finer around the last value)
		 * logarithmic scale (radix=10).
		 */
		INVERSELOGSCALE10;
	}

	/**
	 * The possible types of defining ranges for a {@link RangeParameter}.
	 * 
	 * @author Jan Grau
	 */
	public enum RangeType {
		/**
		 * The parameter is not ranged, i.e. is works like a
		 * {@link SimpleParameter}.
		 */
		NO,
		/**
		 * The values of the parameter are given as a list.
		 */
		LIST,
		/**
		 * The values of the parameter are given as a range of possible values
		 * and a number of steps between the first and the last value.
		 */
		RANGE;
	}

	/**
	 * The {@link SimpleParameter} that is ranged
	 */
	private SimpleParameter rangedParameter;

	/**
	 * The array of parameter values
	 */
	private Object[] values;

	/**
	 * The index of the current value in <code>values</code>
	 */
	private int current;

	/**
	 * The <code>enum</code> that determines the type of range (one of
	 * <code>NO, LIST</code> or <code>RANGE</code>)
	 * 
	 * @see RangeType
	 */
	private RangeType shallBeRanged;

	/**
	 * The <code>enum</code> that determines the type of scale (one of
	 * <code>LOGSCALE, INVERSELOGSCALE, LINSCALE, EXPSCALE,
	 * LOGSCALE2, INVERSELOGSCALE2, LOGSCALE10</code> or
	 * <code>INVERSELOGSCALE10</code>)
	 * 
	 * @see Scale
	 */
	private Scale scale;

	/**
	 * <code>true</code> if the parameters are set
	 */
	private boolean isSet;

	/**
	 * The error message or <code>null</code> if no error occurred
	 */
	private String errorMessage;

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#clone()
	 */
	@Override
	public RangeParameter clone() throws CloneNotSupportedException {
		RangeParameter clone = (RangeParameter) super.clone();
		clone.values = values == null ? null : values.clone();
		clone.rangedParameter = rangedParameter.clone();
		return clone;
	}

	/**
	 * Constructs a {@link RangeParameter} from a {@link SimpleParameter} that
	 * is rangeable.
	 * 
	 * @param par
	 *            the rangeable {@link SimpleParameter}
	 * @throws Exception
	 *             if the {@link SimpleParameter} is not rangeable, i.e. its
	 *             method {@link SimpleParameter#isRangeable()} returns
	 *             <code>false</code>
	 */
	public RangeParameter(SimpleParameter par) throws Exception {
		super( par.getName(), par.getComment(), par.getDatatype() );
		if (!par.isRangeable()) {
			throw new Exception("Parameter must be rangeable");
		}
		this.rangedParameter = par.clone();
		this.shallBeRanged = RangeType.NO;
		this.scale = Scale.LINSCALE;
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Restores a {@link RangeParameter} from its XML representation.
	 * 
	 * @param representation
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the {@link StringBuffer} <code>representation</code> could
	 *             not be parsed
	 */
	public RangeParameter(StringBuffer representation) throws NonParsableException {
		super(representation);
	}

	/**
	 * Returns a {@link CollectionParameter} that allows the user to choose
	 * between different scales.
	 * 
	 * @return the {@link CollectionParameter}
	 * 
	 * @throws ParameterException
	 *             if the {@link CollectionParameter} could, or if one of the
	 *             possible values of the {@link CollectionParameter} was not
	 *             legal
	 */
	public static CollectionParameter getCollectionOfScales()
			throws ParameterException {

		return new EnumParameter(Scale.class,
				"The possible scales for a range of parameter values", true);
	}

	/**
	 * Returns a list of all parameter values as a {@link String} or
	 * <code>null</code> if no parameter values have been set.
	 * 
	 * @return the list of all parameter values or <code>null</code>
	 */
	public String getList() {
		if (values == null) {
			return null;
		} else {
			StringBuffer buf = new StringBuffer();
			for (int i = 0; i < values.length - 1; i++) {
				buf.append(values[i].toString() + " ");
			}
			buf.append(values[values.length - 1].toString());

			return buf.toString();
		}
	}

	/**
	 * Returns the start value of a range of parameter values or
	 * <code>null</code> if no range was specified.
	 * 
	 * @return the start value of a range of parameter values or
	 *         <code>null</code>
	 */
	public Object getStartValue() {
		if (shallBeRanged != RangeType.RANGE) {
			return null;
		} else if (values == null) {
			return rangedParameter.getValue();
		} else {
			return values[0];
		}
	}

	/**
	 * Returns the last value of a range of parameter values or
	 * <code>null</code> if no range was specified.
	 * 
	 * @return the last value of a range of parameter values or
	 *         <code>null</code>
	 */
	public Object getEndValue() {
		if (shallBeRanged != RangeType.RANGE) {
			return null;
		} else if (values == null) {
			return rangedParameter.getValue();
		} else {
			return values[values.length - 1];
		}
	}

	/**
	 * Returns the number of steps of a range of parameter values or
	 * <code>0</code> if no range was specified.
	 * 
	 * @return the number of steps of a range of parameter values or
	 *         <code>0</code>
	 */
	public int getSteps() {
		if (shallBeRanged != RangeType.RANGE) {
			return 0;
		} else {
			if (values == null) {
				return 1;
			} else {
				return values.length - 1;
			}
		}
	}

	/**
	 * Returns a description of the the scale of a range of parameter values.
	 * 
	 * @return the description of the the scale of a range of parameter values
	 * 
	 * @see Scale
	 */
	public String getScale() {
		if (shallBeRanged != RangeType.RANGE) {
			return null;
		} else {
			switch (scale) {
			case LINSCALE:
				return "Linear";
			case EXPSCALE:
				return "Exponential";
			case LOGSCALE:
				return "Logarithmic";
			case INVERSELOGSCALE:
				return "Inverse logarithmic";
			case LOGSCALE2:
				return "Logarithmic (radix=2)";
			case INVERSELOGSCALE2:
				return "Inverse logarithmic (radix=2)";
			case LOGSCALE10:
				return "Logarithmic (radix=10)";
			case INVERSELOGSCALE10:
				return "Inverse logarithmic (radix=10)";
			default:
				return null;
			}
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#isRequired()
	 */
	@Override
	public boolean isRequired() {
		return true;
	}

	private boolean checkSteps(Object steps) {
		if (steps instanceof String) {
			try {
				if (Integer.parseInt((String) steps) > 0) {
					this.errorMessage = null;
					return true;
				} else {
					this.errorMessage = "The number of steps must be at least 1.";
					return false;
				}
			} catch (NumberFormatException e) {
				this.errorMessage = "The number of steps must be an integer.";
				return false;
			}
		} else if (steps instanceof Integer) {
			if (((Integer) steps) > 0) {
				this.errorMessage = null;
				return true;
			} else {
				this.errorMessage = "The number of steps must be at least 1.";
				return false;
			}
		}
		this.errorMessage = "The number of steps is of illegal type "
				+ steps.getClass();
		return false;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#checkValue(java.lang.Object)
	 */
	@Override
	public boolean checkValue(Object value) {
		boolean r = rangedParameter.checkValue(value);
		if (!r && rangedParameter.getErrorMessage() != null) {
			this.errorMessage = rangedParameter.getErrorMessage();
		}
		return r;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#setValue(java.lang.Object)
	 */
	@Override
	public void setValue(Object value) throws IllegalValueException {
		rangedParameter.setValue(value);
	}

	/**
	 * Sets a list of values from a {@link String} containing a space separated
	 * list of values. After this method has been called,
	 * {@link #shallBeRanged()} returns {@link RangeType#LIST}.
	 * 
	 * @param values
	 *            the space separated list of values
	 * 
	 * @throws IllegalValueException
	 *             if at least one of the values in the list is not suitable for
	 *             this {@link RangeParameter}
	 */
	public void setValues(String values) throws IllegalValueException {
		shallBeRanged = RangeType.LIST;
		StringTokenizer tok = new StringTokenizer(values);
		this.values = new Object[tok.countTokens()];
		for (int i = 0; i < this.values.length; i++) {
			String temp = tok.nextToken();
			try {
				switch (getDatatype()) {
				case BYTE:
					this.values[i] = new Byte(temp);
					break;
				case SHORT:
					this.values[i] = new Short(temp);
					break;
				case INT:
					this.values[i] = new Integer(temp);
					break;
				case LONG:
					this.values[i] = new Long(temp);
					break;
				case FLOAT:
					this.values[i] = new Float(temp);
					break;
				case DOUBLE:
					this.values[i] = new Double(temp);
					break;

				}
			} catch (NumberFormatException e) {
				this.values = null;
				errorMessage = e.getMessage();
				throw new IllegalValueException(e.getMessage());
			}
			if (!rangedParameter.checkValue(this.values[i])) {
				errorMessage = "Value " + temp + " not allowed! "
						+ rangedParameter.getErrorMessage();
				throw new IllegalValueException("Value " + temp
						+ " not allowed! " + rangedParameter.getErrorMessage());
			}
		}
		current = 0;
		setValue(this.values[0]);
		errorMessage = "";
	}

	/**
	 * Sets the values of this {@link RangeParameter} as a range of values,
	 * specified by a start value, a last value, a number of steps between these
	 * values (without the last value) and a scale in that the values between
	 * the first and the last value are chosen. This range must be one of
	 * <code>LOGSCALE, INVERSELOGSCALE, LINSCALE, EXPSCALE,
	 * LOGSCALE2, INVERSELOGSCALE2, LOGSCALE10</code> or
	 * <code>INVERSELOGSCALE10</code>.
	 * 
	 * @param startValue
	 *            the first value of the range
	 * @param steps
	 *            the number of steps between <code>startValue</code> and
	 *            <code>endValue</code> (exclusive)
	 * @param endValue
	 *            the last value of the range
	 * @param scale
	 *            the {@link Scale}
	 * @throws IllegalValueException
	 *             if one of the values is not suitable for this
	 *             {@link RangeParameter}, e.g. is out of a specified range of
	 *             allowed values
	 * 
	 * @see Scale
	 */
	public void setValues(Object startValue, int steps, Object endValue,
			Scale scale) throws IllegalValueException {
		if (scale == Scale.LINSCALE || scale == Scale.EXPSCALE) {
			values = check(startValue, steps, endValue);

			DataType d = getDatatype();
			double stepSize, start = ((Number) values[0]).doubleValue(), end = ((Number) values[steps])
					.doubleValue();

			if (scale == Scale.LINSCALE) {
				stepSize = (double) (end - start) / (double) steps;
				for (int i = 1; i < steps; i++) {
					switch (d) {
					case BYTE:
						values[i] = new Byte((byte) (start + (i * stepSize)));
						break;
					case SHORT:
						values[i] = new Short((short) (start + (i * stepSize)));
						break;
					case INT:
						values[i] = new Integer((int) (start + (i * stepSize)));
						break;
					case LONG:
						values[i] = new Long((long) (start + (i * stepSize)));
						break;
					case FLOAT:
						values[i] = new Float((float) (start + (i * stepSize)));
						break;
					case DOUBLE:
						values[i] = new Double(
								(double) (start + (i * stepSize)));
						break;
					}
				}
			} else {
				// TODO check ... ???
				stepSize = (StrictMath.exp(end) - StrictMath.exp(start))
						/ (double) steps;
				for (int i = 1; i < steps; i++) {
					switch (d) {
					case BYTE:
						values[i] = new Byte((byte) (StrictMath.log(StrictMath
								.exp(start)
								+ (i * stepSize))));
						break;
					case SHORT:
						values[i] = new Short((short) (StrictMath
								.log(StrictMath.exp(start) + (i * stepSize))));
						break;
					case INT:
						values[i] = new Integer((int) (StrictMath
								.log(StrictMath.exp(start) + (i * stepSize))));
						break;
					case LONG:
						values[i] = new Long((long) (StrictMath.log(StrictMath
								.exp(start)
								+ (i * stepSize))));
						break;
					case FLOAT:
						values[i] = new Float((float) (StrictMath
								.log(StrictMath.exp(start) + (i * stepSize))));
						break;
					case DOUBLE:
						values[i] = new Double((double) (StrictMath
								.log(StrictMath.exp(start) + (i * stepSize))));
						break;
					}
				}
			}

			current = 0;
			setValue(values[0]);
		} else {
			boolean log = true;
			double radix = 0;
			switch (scale) {
			case INVERSELOGSCALE:
			case INVERSELOGSCALE2:
			case INVERSELOGSCALE10:
				log = false;
			case LOGSCALE:
			case LOGSCALE2:
			case LOGSCALE10:
				switch (scale) {
				case INVERSELOGSCALE:
				case LOGSCALE:
					radix = Math.E;
					break;
				case INVERSELOGSCALE2:
				case LOGSCALE2:
					radix = 2;
					break;
				case INVERSELOGSCALE10:
				case LOGSCALE10:
					radix = 10;
					break;
				}
				setValuesInLogScale(log, radix, startValue, steps, endValue);
				break;
			default:
				errorMessage = "Range type not allowed!";
				values = null;
				throw new IllegalValueException("Range not allowed!");
			}
		}

		shallBeRanged = RangeType.RANGE;
		this.scale = scale;

		errorMessage = "";
	}

	/**
	 * Checks some constraints and returns an array that can be used for storing
	 * all values containing start and end value of the range.
	 * 
	 * @param startValue
	 *            the first value of the range
	 * @param steps
	 *            the number of steps between <code>startValue</code> and
	 *            <code>endValue</code> (exclusive)
	 * @param endValue
	 *            the last value of the range
	 * 
	 * @return an array that can be used for storing all values, containing
	 *         already the start and end value
	 * 
	 * @throws IllegalValueException
	 *             if some constraints are not fulfilled
	 */
	private Object[] check(Object startValue, int steps, Object endValue)
			throws IllegalValueException {
		if (steps < 1) {
			errorMessage = "The number of steps has to be at least 1.";
			throw new IllegalValueException(
					"The number of steps has to be at least 1.");
		}
		if (!rangedParameter.checkValue(startValue)) {
			errorMessage = "Start value out of range!";
			throw new IllegalValueException("Start value out of range!");
		}
		if (!rangedParameter.checkValue(endValue)) {
			errorMessage = "End value out of range!";
			throw new IllegalValueException("End value out of range!");
		}
		Object[] val = new Object[steps + 1];
		rangedParameter.setValue(startValue);
		val[0] = rangedParameter.getValue();
		rangedParameter.setValue(endValue);
		val[steps] = rangedParameter.getValue();
		if (((Number) val[0]).doubleValue() > ((Number) val[steps])
				.doubleValue()) {
			errorMessage = "Start value has to be less than end value!";
			throw new IllegalValueException(
					"Start value has to be less than end value!");
		}
		return val;
	}

	/**
	 * This method enables you to set a list of values in an easy manner. The
	 * values of the list are from the interval
	 * <code>[startValue,endValue]</code> in (inverse) logarithmic scale.
	 * 
	 * @param log
	 *            if <code>true</code> a logarithmic scale is used, otherwise
	 *            the inverse logarithmic scale
	 * @param radix
	 *            the radix for the scale
	 * @param startValue
	 *            the first value of the range
	 * @param steps
	 *            the number of steps between <code>startValue</code> and
	 *            <code>endValue</code> (exclusive)
	 * @param endValue
	 *            the last value of the range
	 * 
	 * @throws IllegalValueException
	 *             if some constraints are not fulfilled
	 */
	public void setValuesInLogScale(boolean log, double radix,
			Object startValue, int steps, Object endValue)
			throws IllegalValueException {
		DataType d = getDatatype();
		values = check(startValue, steps, endValue);
		double start = ((Number) values[0]).doubleValue(), end = ((Number) values[steps])
				.doubleValue(), stepSize = end - start;
		if (log) {
			for (int i = 1; i < steps; i++) {
				switch (d) {
				case BYTE:
					values[i] = new Byte((byte) (start + Math.pow(radix, i
							- steps)
							* stepSize));
					break;
				case SHORT:
					values[i] = new Short((short) (start + Math.pow(radix, i
							- steps)
							* stepSize));
					break;
				case INT:
					values[i] = new Integer((int) (start + Math.pow(radix, i
							- steps)
							* stepSize));
					break;
				case LONG:
					values[i] = new Long((long) (start + Math.pow(radix, i
							- steps)
							* stepSize));
					break;
				case FLOAT:
					values[i] = new Float((float) (start + Math.pow(radix, i
							- steps)
							* stepSize));
					break;
				case DOUBLE:
					values[i] = new Double((double) (start + Math.pow(radix, i
							- steps)
							* stepSize));
					break;
				}
			}
		} else {
			for (int i = 1; i < steps; i++) {
				switch (d) {
				case BYTE:
					values[steps - i] = new Byte((byte) (end - Math.pow(radix,
							i - steps)
							* stepSize));
					break;
				case SHORT:
					values[steps - i] = new Short((short) (end - Math.pow(
							radix, i - steps)
							* stepSize));
					break;
				case INT:
					values[steps - i] = new Integer((int) (end - Math.pow(
							radix, i - steps)
							* stepSize));
					break;
				case LONG:
					values[steps - i] = new Long((long) (end - Math.pow(radix,
							i - steps)
							* stepSize));
					break;
				case FLOAT:
					values[steps - i] = new Float((float) (end - Math.pow(
							radix, i - steps)
							* stepSize));
					break;
				case DOUBLE:
					values[steps - i] = new Double((double) (end - Math.pow(
							radix, i - steps)
							* stepSize));
					break;
				}
			}
		}
		current = 0;
		setValue(values[0]);
		shallBeRanged = RangeType.LIST;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#getValue()
	 */
	@Override
	public Object getValue() {
		return rangedParameter.getValue();
	}

	/**
	 * Returns <code>true</code> if the next element still exists and can be
	 * fetched using {@link #getValue()}, <code>false</code> otherwise.
	 * 
	 * @return <code>true</code> if the next element exists and can be fetched,
	 *         <code>false</code> otherwise
	 * 
	 * @throws IllegalValueException
	 *             if the next value to be fetched is not valid
	 * 
	 * @see RangeParameter#setValue(Object)
	 */
	public boolean next() throws IllegalValueException {
		if (values != null && current + 1 < values.length) {
			current++;
			setValue(values[current]);
			return true;
		} else {
			return false;
		}
	}

	/**
	 * Returns the number of values in a list or range of parameter values.
	 * 
	 * @return the number of values
	 */
	public int getNumberOfValues() {
		if (values != null) {
			return values.length;
		} else {
			return 1;
		}
	}

	/**
	 * Returns the number of calls of {@link #next()} that can be done before
	 * obtaining <code>false</code>.
	 * 
	 * @param afterIndex
	 *            the index after which {@link #next()} shall be called
	 * 
	 * @return the number of possible calls of {@link #next()} after
	 *         <code>afterIndex</code>
	 */
	public int getNumberOfNexts(int afterIndex) {
		if (values != null) {
			return values.length - afterIndex - 1;
		} else {
			return 0;
		}
	}

	/**
	 * Returns one of <code>LIST, RANGE</code> or <code>NO</code> depending on
	 * the input used to specify this {@link RangeParameter}.
	 * 
	 * @return the type of range specification
	 * 
	 * @see RangeType
	 */
	public RangeType shallBeRanged() {
		return shallBeRanged;
	}

	/**
	 * Sets the type of this {@link RangeParameter} to one of
	 * <code>LIST, RANGE</code> or <code>NO</code>.
	 * 
	 * @param shallBeRanged
	 *            the {@link RangeType} of this parameter
	 * 
	 * @throws Exception
	 *             if <code>shallBeRanged</code> is none of the allowed values
	 * 
	 * @see RangeType
	 */
	public void setShallBeRanged(RangeType shallBeRanged) throws Exception {
		if (shallBeRanged == RangeType.NO || shallBeRanged == RangeType.LIST
				|| shallBeRanged == RangeType.RANGE) {
			this.shallBeRanged = shallBeRanged;
			if (shallBeRanged == RangeType.NO) {
				this.values = null;
			}
		} else {
			throw new Exception("None of the range parameters!");
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#hasDefaultOrIsSet()
	 */
	@Override
	public boolean hasDefaultOrIsSet() {
		if (shallBeRanged == RangeType.NO) {
			return rangedParameter.hasDefaultOrIsSet();
		} else {
			return (values != null);
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#isSet()
	 */
	@Override
	public boolean isSet() {
		if (shallBeRanged == RangeType.NO) {
			return rangedParameter.isSet();
		} else {
			return (values != null);
		}
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
	 * @see de.jstacs.parameters.Parameter#getErrorMessage()
	 */
	@Override
	public String getErrorMessage() {
		if (errorMessage == null) {
			return rangedParameter.getErrorMessage();
		} else {
			return errorMessage;
		}
	}
	
	/*
	 * (non-Javadoc)
	 * @see de.jstacs.AnnotatedEntity#getXMLTag()
	 */
	@Override
	public String getXMLTag() {
		return "rangeParameter";
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.parameters.Parameter#appendFurtherInfos(java.lang.StringBuffer)
	 */
	@Override
	protected void appendFurtherInfos( StringBuffer buf ) {
		super.appendFurtherInfos( buf );
		
		XMLParser.appendObjectWithTags(buf, rangedParameter, "rangedParameter");
		if (values != null) {
			StringBuffer buf2 = new StringBuffer();
			XMLParser.appendObjectWithTags(buf2, values, "vals");
			XMLParser.appendObjectWithTags(buf, buf2.toString(), "values");
		} else {
			XMLParser.appendObjectWithTags(buf, null, "values");
		}
		XMLParser.appendObjectWithTags(buf, current, "current");
		XMLParser.appendObjectWithTags(buf, isSet, "isSet");
		XMLParser.appendObjectWithTags(buf, shallBeRanged, "shallBeRanged");
		XMLParser.appendObjectWithTags(buf, errorMessage, "errorMessage");
		XMLParser.appendObjectWithTags(buf, scale, "scale");

	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#fromXML(java.lang.StringBuffer)
	 */
	@Override
	protected void extractFurtherInfos( StringBuffer buf ) throws NonParsableException {
		super.extractFurtherInfos( buf );
		rangedParameter = XMLParser.extractObjectForTags(buf,"rangedParameter", SimpleParameter.class);
		StringBuffer buf2 = XMLParser.extractForTag(buf, "values");
		if (buf2.toString().equals("null")) {
			values = null;
		} else {
			values = (Object[]) XMLParser.extractObjectForTags(buf2, "vals");
		}
		current = XMLParser.extractObjectForTags(buf, "current", int.class );
		isSet = XMLParser.extractObjectForTags(buf, "isSet", boolean.class );
		shallBeRanged = XMLParser.extractObjectForTags(buf, "shallBeRanged", RangeType.class );
		errorMessage = XMLParser.parseString( XMLParser.extractObjectForTags(buf, "errorMessage", String.class ) );
		scale = XMLParser.extractObjectForTags(buf, "scale", Scale.class );
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#reset()
	 */
	@Override
	public void reset() {
		rangedParameter.reset();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.RangeIterator#resetToFirst()
	 */
	public void resetToFirst() {
		current = 0;
		try {
			setValue(values[current]);
		} catch (Exception e) {
			rangedParameter.reset();
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.Parameter#setDefault(java.lang.Object)
	 */
	@Override
	public void setDefault(Object defaultValue) throws Exception {
		rangedParameter.setDefault(defaultValue);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.RangeIterator#valuesToString()
	 */
	public String valuesToString() {
		if (shallBeRanged == RangeType.LIST) {
			return "[List: " + getList() + "]";
		} else if (shallBeRanged == RangeType.RANGE) {
			return "[" + getScale() + ": " + getStartValue() + " .. "
					+ getSteps() + " .. " + getEndValue() + "]";
		} else {
			return "[" + getValue().toString() + "]";
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.parameters.RangeIterator#isRanged()
	 */
	public boolean isRanged() {
		return shallBeRanged != RangeType.NO;
	}


	@Override
	public void toGalaxy( String namePrefix, String configPrefix, int depth, StringBuffer descBuffer, StringBuffer configBuffer ) throws Exception {
		
		namePrefix = namePrefix+"_"+GalaxyAdaptor.getLegalName( getName() );
		configPrefix = configPrefix+namePrefix+".";
		
		StringBuffer buf = new StringBuffer();
		
		EnumParameter temp = new EnumParameter( RangeType.class, "Define a single value (NO), a list of values (LIST) or a range of values (RANGE) for "+rangedParameter.getName(), true );
		temp.toGalaxy( namePrefix, configPrefix, depth+1, buf, configBuffer );
		
		
		StringBuffer tmp = new StringBuffer();
		configBuffer.append( "#if $"+configPrefix+namePrefix+"_RangeType == \"NO\"\n" );
		rangedParameter.toGalaxy( namePrefix+"_NO", configPrefix, depth+1, tmp, configBuffer );
		XMLParser.addTagsAndAttributes( tmp, "when", "value=\"NO\"");
		buf.append( tmp );
		
		tmp = new StringBuffer();
		configBuffer.append( "#elif $"+configPrefix+namePrefix+"_RangeType == \"LIST\"\n" );
		(new SimpleParameter( DataType.STRING, rangedParameter.getName(), "Please enter a comma-separated list of values. "+rangedParameter.getComment(), true )).toGalaxy( namePrefix+"_LIST", configPrefix, depth+1, tmp, configBuffer );
		XMLParser.addTagsAndAttributes( tmp, "when", "value=\"LIST\"");
		buf.append( tmp );
		
		tmp = new StringBuffer();
		configBuffer.append( "#elif $"+configPrefix+namePrefix+"_RangeType == \"RANGE\"\n" );
		(new SimpleParameter( rangedParameter.getDatatype(), rangedParameter.getName()+": start value", "Please enter the start value. "+rangedParameter.getComment(), true )).toGalaxy( namePrefix+"_RANGE", configPrefix, depth+1, tmp, configBuffer);
		(new SimpleParameter( rangedParameter.getDatatype(), rangedParameter.getName()+": end value", "Please enter the end value. "+rangedParameter.getComment(), true )).toGalaxy( namePrefix+"_RANGE", configPrefix, depth+1, tmp, configBuffer );
		(new SimpleParameter( rangedParameter.getDatatype(), rangedParameter.getName()+": steps", "Please enter the number of steps between start and end value. "+rangedParameter.getComment(), true )).toGalaxy( namePrefix+"_RANGE", configPrefix, depth+1, tmp, configBuffer );
		temp = new EnumParameter( Scale.class, "Select a scaling between start and end value.", true );
		temp.toGalaxy( namePrefix+"_RANGE", configPrefix, depth+1, tmp, configBuffer );
		XMLParser.addTagsAndAttributes( tmp, "when", "value=\"RANGE\"" );
		buf.append( tmp );
		configBuffer.append( "#end if\n" );
		
		XMLParser.addTagsAndAttributes( buf, "conditional", "name=\""+namePrefix+"\"" );
		
		descBuffer.append( buf );		
	}

	@Override
	public void fromGalaxy( String namePrefix, StringBuffer command ) throws Exception {
		namePrefix = namePrefix+"_"+GalaxyAdaptor.getLegalName( getName() );
		String type = XMLParser.extractForTag( command, namePrefix+"_RangeType" ).toString();
		this.setShallBeRanged( RangeType.valueOf( type ) );
		if(this.shallBeRanged == RangeType.NO){
			this.setValue( XMLParser.extractForTag( command, namePrefix+"_NO_"+rangedParameter.getName().replaceAll( "[\\s:-]+", "_" ) ).toString() );
		}else if(this.shallBeRanged == RangeType.LIST){
			this.setValues( XMLParser.extractForTag( command, namePrefix+"_LIST_"+rangedParameter.getName().replaceAll( "[\\s:-]+", "_" ) ).toString() );
		}else{
			this.setValues( XMLParser.extractForTag( command, namePrefix+"_RANGE_"+rangedParameter.getName().replaceAll( "[\\s:-]+", "_" )+"_start_value" ).toString(),
					XMLParser.extractObjectForTags( command, namePrefix+"_RANGE_"+rangedParameter.getName().replaceAll( "[\\s:-]+", "_" )+"_steps", int.class ),
					XMLParser.extractForTag( command, namePrefix+"_RANGE_"+rangedParameter.getName().replaceAll( "[\\s:-]+", "_" )+"_end_value" ).toString(),
					Scale.valueOf( XMLParser.extractForTag( command, namePrefix+"_RANGE_Scale" ).toString() ) );
		}
		
		
	}
	
}