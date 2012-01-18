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

package de.jstacs.sequenceScores.statisticalModels.trainable;

import java.io.IOException;
import java.util.LinkedList;

import de.jstacs.NotTrainedException;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.io.ArrayHandler;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.results.NumericalResult;
import de.jstacs.results.NumericalResultSet;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.results.StorableResult;

/**
 * This class is for modelling sequences by modelling the different positions of
 * the each sequence by different models. For instance one can use this class to
 * model a sequence of subsequences, where each subsequence is modeled by a
 * different model.
 * 
 * @author Jens Keilwagen
 */
public class CompositeTrainSM extends AbstractTrainSM {
	private static final long serialVersionUID = 1263707296720984521L;

	private static final String XML_TAG = "CompositeTrainSM";

	/**
	 * The models for the components
	 */
	protected TrainableStatisticalModel[] models;

	/**
	 * The start indices.
	 */
	protected int[][] starts;

	/**
	 * The length for each component.
	 */
	protected int[][] lengths;

	/**
	 * Creates a new {@link CompositeTrainSM}.
	 * 
	 * @param alphabets
	 *            the alphabets used in the models
	 * @param assignment
	 *            an array assigning each position to a model
	 * @param models
	 *            the single models
	 * 
	 * @throws IllegalArgumentException
	 *             if something is wrong with the assignment of the models
	 * @throws WrongAlphabetException
	 *             if (parts of) the alphabet does not match with the alphabets
	 *             of the models
	 * @throws CloneNotSupportedException
	 *             if the models could not be cloned
	 */
	public CompositeTrainSM(AlphabetContainer alphabets, int[] assignment,
			TrainableStatisticalModel... models) throws WrongAlphabetException,
			CloneNotSupportedException {
		super(alphabets, assignment.length);
		starts = new int[models.length][];
		lengths = new int[models.length][];
		int i = 0, start = -1;
		int[] len = new int[models.length], parts = new int[models.length];
		for (; i < assignment.length; i++) {
			if (assignment[i] != start) {
				parts[assignment[i]]++;
				start = assignment[i];
			}
		}
		for (i = 0; i < models.length; i++) {
			starts[i] = new int[parts[i]];
			lengths[i] = new int[parts[i]];
			parts[i] = 0;
		}
		len[assignment[0]] = 1;
		for (start = 0, i = 1; i < assignment.length; i++) {
			len[assignment[i]]++;
			if (assignment[i] != assignment[start]) {
				// another model shall be used

				// set last subsequence
				starts[assignment[start]][parts[assignment[start]]] = start;
				lengths[assignment[start]][parts[assignment[start]]] = i
						- start;
				parts[assignment[start]]++;
				// set new model
				start = i;
			}
		}
		starts[assignment[start]][parts[assignment[start]]] = start;
		lengths[assignment[start]][parts[assignment[start]]] = i - start;

		// check alphabets and lengths
		for (i = 0; i < models.length; i++) {
			if (!alphabets.getCompositeContainer(starts[i], lengths[i])
					.checkConsistency(models[i].getAlphabetContainer())) {
				throw new WrongAlphabetException("The " + i
						+ "-th model has not the correct alphabet.");
			}
			if (models[i].getLength() != 0 && models[i].getLength() != len[i]) {
				throw new IllegalArgumentException("The " + i
						+ "-th model can not model sequence of length "
						+ len[i] + ", model length is " + models[i].getLength()
						+ ".");
			}
		}
		this.models = ArrayHandler.clone(models);
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Creates a new {@link CompositeTrainSM} out of a {@link StringBuffer}.
	 * 
	 * @param stringBuff
	 *            the {@link StringBuffer} to be parsed
	 * 
	 * @throws NonParsableException
	 *             if the {@link StringBuffer} is not parsable
	 */
	public CompositeTrainSM(StringBuffer stringBuff) throws NonParsableException {
		super(stringBuff);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.trainableStatisticalModels.AbstractTrainSM#clone()
	 */
	@Override
	public CompositeTrainSM clone() throws CloneNotSupportedException {
		CompositeTrainSM c = (CompositeTrainSM) super.clone();
		c.starts = new int[models.length][];
		c.lengths = new int[models.length][];
		for (int i = 0; i < models.length; i++) {
			c.starts[i] = starts[i].clone();
			c.lengths[i] = lengths[i].clone();
		}
		c.models = ArrayHandler.clone(models);
		return c;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.trainableStatisticalModels.AbstractTrainSM#getCharacteristics()
	 */
	@Override
	public ResultSet getCharacteristics() throws Exception {
		LinkedList<Result> infos = new LinkedList<Result>();
		ResultSet part;
		int i = 0, j;
		for (; i < models.length; i++) {
			part = models[i].getCharacteristics();
			if (part != null) {
				infos
						.add(new NumericalResult("model number",
								"type of model "
										+ models[i].getClass().getSimpleName(),
								new Integer(i)));
				for (j = 0; j < part.getNumberOfResults(); j++) {
					infos.add(part.getResultAt(j));
				}
			}
		}
		infos.add(new StorableResult("model",
				"the xml representation of the model", this));
		return new ResultSet(infos);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.trainableStatisticalModels.TrainableStatisticalModel#getInstanceName()
	 */
	public String getInstanceName() {
		StringBuffer erg = new StringBuffer("composite model(");
		erg.append(models[0].getInstanceName());
		for (int i = 1; i < models.length; i++) {
			erg.append(", ");
			erg.append(models[i].getInstanceName());
		}
		erg.append(")");
		return erg.toString();
	}

	/**
	 * This method returns the length of the models in the
	 * {@link CompositeTrainSM}.
	 * 
	 * @return the length of the models as an array
	 */
	public int[] getLengthOfModels() {
		int[] erg = new int[models.length];
		for (int i = 0; i < models.length; i++) {
			erg[i] = models[i].getLength();
		}
		return erg;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.trainableStatisticalModels.AbstractTrainSM#getMaximalMarkovOrder()
	 */
	@Override
	public byte getMaximalMarkovOrder() throws UnsupportedOperationException {
		byte max = 0, c;
		for (int i = 0; i < models.length; i++) {
			c = models[i].getMaximalMarkovOrder();
			if (max < c) {
				max = c;
			}
		}
		return max;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.trainableStatisticalModels.TrainableStatisticalModel#getNumericalCharacteristics()
	 */
	public NumericalResultSet getNumericalCharacteristics() throws Exception {
		LinkedList<NumericalResult> infos = new LinkedList<NumericalResult>();
		NumericalResultSet part;
		int i = 0, j;
		for (; i < models.length; i++) {
			part = models[i].getNumericalCharacteristics();
			if (part != null && part.getNumberOfResults() > 0) {
				infos
						.add(new NumericalResult("model number",
								"type of model "
										+ models[i].getClass().getSimpleName(),
								new Integer(i)));
				for (j = 0; j < part.getNumberOfResults(); j++) {
					infos.add(part.getResultAt(j));
				}
			}
		}
		return new NumericalResultSet(infos);
	}

	/**
	 * Returns the a deep copy of the models.
	 * 
	 * @return an array of {@link AbstractTrainSM}s
	 * 
	 * @throws CloneNotSupportedException
	 *             if at least one of the models could not be cloned
	 */
	public TrainableStatisticalModel[] getModels() throws CloneNotSupportedException {
		return ArrayHandler.clone(models);
	}

	/**
	 * This method returns the number of models in the {@link CompositeTrainSM}.
	 * 
	 * @return the number of models
	 */
	public int getNumberOfModels() {
		return models.length;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.trainableStatisticalModels.TrainableStatisticalModel#getLogPriorTerm()
	 */
	public double getLogPriorTerm() throws Exception {
		double sum = 0;
		for (int i = 0; i < models.length; i++) {
			sum += models[i].getLogPriorTerm();
		}
		return sum;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * de.jstacs.trainableStatisticalModels.AbstractTrainSM#getLogProbFor(de.jstacs.data.Sequence,
	 * int, int)
	 */
	@Override
	public double getLogProbFor(Sequence sequence, int startpos, int endpos)
			throws NotTrainedException, Exception {
		if (endpos - startpos + 1 != length) {
			throw new IllegalArgumentException("This sequence has not length "
					+ length + ".");
		}
		double erg = 0;
		for (int i = 0; i < models.length; i++) {
			if (lengths[i].length == 1) {
				erg += models[i].getLogProbFor(sequence, startpos
						+ starts[i][0], startpos + starts[i][0] + lengths[i][0]
						- 1);
			} else {
				erg += models[i].getLogProbFor(sequence.getCompositeSequence(
						alphabets, getStartsFor(i, startpos), lengths[i]));
			}
		}
		return erg;
	}

	private int[] getStartsFor(int index, int offset) {
		int[] arrays = new int[starts[index].length];
		for (int i = 0; i < arrays.length; i++) {
			arrays[i] = offset + starts[index][i];
		}
		return arrays;
	}

	/*
	 * (non-Javadoc)
	 * @see de.jstacs.SequenceScore#isInitialized()
	 */
	public boolean isInitialized() {
		boolean erg = true;
		int i = 0;
		while (i < models.length && erg) {
			erg &= models[i++].isInitialized();
		}
		return erg;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.trainableStatisticalModels.AbstractTrainSM#fromXML(java.lang.StringBuffer)
	 */
	@Override
	public void fromXML(StringBuffer representation)
			throws NonParsableException {
		StringBuffer xml = XMLParser.extractForTag(representation, XML_TAG);
		alphabets = new AlphabetContainer(xml);
		length = XMLParser.extractObjectForTags(xml, "sequencelength", int.class );
		lengths = XMLParser.extractObjectForTags(xml, "lengths", int[][].class );
		starts = XMLParser.extractObjectForTags(xml, "starts", int[][].class );
		models = XMLParser.extractObjectForTags( xml, "models", AbstractTrainSM[].class );
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.Storable#toXML()
	 */
	public StringBuffer toXML() {
		StringBuffer xml = alphabets.toXML();
		XMLParser.appendObjectWithTags(xml, length, "sequencelength");
		XMLParser.appendObjectWithTags(xml, lengths, "lengths");
		XMLParser.appendObjectWithTags(xml, starts, "starts");
		XMLParser.appendObjectWithTags(xml, models, "models");
		XMLParser.addTags(xml, XML_TAG);
		return xml;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.trainableStatisticalModels.TrainableStatisticalModel#train(de.jstacs.data.DataSet, double[])
	 */
	public void train(DataSet data, double[] weights) throws Exception {
		if (data.getElementLength() != length) {
			throw new IOException(
					"The given data has not correct sequence length.");
		}
		for (int i = 0; i < models.length; i++) {
			if (weights == null) {
				models[i].train(data.getCompositeDataSet(starts[i], lengths[i]));
			} else {
				models[i].train(data.getCompositeDataSet(starts[i], lengths[i]), weights);
			}
		}
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		String erg = "";
		for (int i = 0; i < models.length; i++) {
			erg += "model: " + i + "\n";
			erg += models[i] + "\n\n";
		}
		return erg;
	}
}
