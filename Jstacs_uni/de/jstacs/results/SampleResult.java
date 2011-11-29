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

package de.jstacs.results;

import de.jstacs.DataType;
import de.jstacs.NonParsableException;
import de.jstacs.Storable;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.Sample;
import de.jstacs.data.Sequence;
import de.jstacs.data.Sample.ElementEnumerator;
import de.jstacs.data.sequences.annotation.SequenceAnnotation;
import de.jstacs.data.sequences.annotation.SequenceAnnotationParser;
import de.jstacs.io.ArrayHandler;
import de.jstacs.io.StringExtractor;
import de.jstacs.io.XMLParser;

/**
 * {@link Result} that contains a {@link Sample}. This {@link Sample} e.g. may
 * have been created by {@link de.jstacs.models.Model#emitSample(int, int...)},
 * or maybe a {@link Sample} that has been annotated in a classification.
 * 
 * @author Jan Grau
 */
public class SampleResult extends Result {

	private Sample data;
	private SequenceAnnotationParser parser;

	/**
	 * Creates a new {@link SampleResult} from a {@link Sample} with the
	 * annotation <code>name</code> and <code>comment</code>.
	 * 
	 * @param name
	 *            the name of the {@link Result}
	 * @param comment
	 *            the comment on the {@link Result}
	 * @param data
	 *            the {@link Sample} that is the result of some computation
	 */
	public SampleResult(String name, String comment, Sample data) {
		super(name, comment, DataType.SAMPLE);
		this.data = data;
	}
	
	/**
	 * Creates a new {@link SampleResult} from a {@link Sample} with the
	 * annotation <code>name</code> and <code>comment</code>.
	 * 
	 * @param name
	 *            the name of the {@link Result}
	 * @param comment
	 *            the comment on the {@link Result}
	 * @param data
	 *            the {@link Sample} that is the result of some computation
	 * @param parser 
	 * 			  a {@link SequenceAnnotationParser} that can be used to store and parse the annotations of <code>data</code>
	 */
	public SampleResult(String name, String comment, Sample data, SequenceAnnotationParser parser) {
		this(name, comment, data);
		this.parser = parser;
	}

	/**
	 * The standard constructor for the interface {@link de.jstacs.Storable}.
	 * Recreates a {@link SampleResult} from its XML representation as returned
	 * by {@link #toXML()}.
	 * 
	 * @param source
	 *            the XML representation as {@link StringBuffer}
	 * 
	 * @throws NonParsableException
	 *             if the XML code could not be parsed
	 */
	public SampleResult(StringBuffer source) throws NonParsableException {
		super(source);
	}
	
	

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.results.Result#fromXML(java.lang.StringBuffer)
	 */
	@Override
	protected void fromXML(StringBuffer rep) throws NonParsableException {
		rep = XMLParser.extractForTag(rep, "sampleResult");
		AlphabetContainer cont = XMLParser.extractObjectForTags(rep, "alphabet", AlphabetContainer.class );// TODO XMLP14CONV This and (possibly) the following lines have been converted automatically
		String sampleAnn = XMLParser.extractObjectForTags(rep, "sampleAnnotation", String.class );
		String seqs = XMLParser.extractObjectForTags(rep, "data", String.class );// TODO XMLP14CONV This and (possibly) the following lines have been converted automatically
		StringExtractor ex = new StringExtractor(seqs, 100, "");
		try {
			Sample data = new Sample(cont, ex);
			Sequence[] seq = data.getAllElements();

			Storable[][] annotation = XMLParser.extractObjectForTags( rep, "annotation", Storable[][].class);
			for (int i = 0; i < seq.length; i++) {
				if (annotation[i].length > 0) {
					seq[i] = seq[i].annotate(false,
							ArrayHandler.cast(SequenceAnnotation.class, annotation[i]));
				}
			}

			this.data = new Sample(sampleAnn, seq);
		} catch (Exception e) {
			NonParsableException np = new NonParsableException(e.getCause()
					.getMessage());
			np.setStackTrace(e.getStackTrace());
			throw np;
		}

	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.results.Result#getResult()
	 */
	@Override
	public Sample getResult() {
		return data;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.Storable#toXML()
	 */
	public StringBuffer toXML() {
		StringBuffer buf = new StringBuffer();
		AlphabetContainer cont = data.getAlphabetContainer();
		XMLParser.appendObjectWithTags(buf, cont, "alphabet");
		StringBuffer seqs = new StringBuffer();
		Storable[][] annotation = new Storable[data.getNumberOfElements()][];
		ElementEnumerator it = new ElementEnumerator(data);
		int i = 0;
		while (it.hasMoreElements()) {
			Sequence seq = it.nextElement();
			buf.append(seq.toString());
			buf.append("\n");
			SequenceAnnotation[] ann = seq.getAnnotation();
			if (ann == null) {
				annotation[i] = new Storable[0];
			} else {
				annotation[i] = ann;
			}
			i++;
		}
		XMLParser.appendObjectWithTags(buf, data.getAnnotation(),
				"sampleAnnotation");
		XMLParser.appendObjectWithTags(buf, seqs.toString(), "data");
		XMLParser.appendObjectWithTags(buf, annotation, "annotation");
		XMLParser.addTags(buf, "sampleResult");
		return buf;
	}

	/**
	 * Returns the {@link SequenceAnnotationParser} that can be used to
	 * write this {@link SampleResult} including annotations on the contained {@link Sequence}s
	 * to a file.
	 * @return the {@link SequenceAnnotationParser}
	 */
	public SequenceAnnotationParser getParser() {
		return parser;
	}

	/**
	 * Sets the {@link SequenceAnnotationParser} that can be used to
	 * write this {@link SampleResult} including annotations on the contained {@link Sequence}s
	 * to a file
	 * @param parser the new {@link SequenceAnnotationParser}
	 */
	public void setParser( SequenceAnnotationParser parser ) {
		this.parser = parser;
	}

}
