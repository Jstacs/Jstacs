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

package de.jstacs.motifDiscovery;

import de.jstacs.data.sequences.Sequence;
import de.jstacs.motifDiscovery.MotifDiscoverer.KindOfProfile;
import de.jstacs.results.ImageResult;
import de.jstacs.utils.REnvironment;

/**
 * This class contains static methods for the {@link MotifDiscoverer}.
 * 
 * @author Jens Keilwagen
 * 
 * @see de.jstacs.motifDiscovery.MotifDiscoverer
 */
public class MotifDiscovererToolBox {

	private final static String PLOT = "plot( 0:(length(profile)-1),profile, \"l\", ylab=\"score\", xlab=\"position\", ylim=lim );\n";

	private final static String ANNOTATE = "abline( h = threshold, lty=2, col=2 );\n"
			+ "for( i in 1:length(profile)){\n"
			+ "if( profile[i] >= threshold ){ h = paste( i-1, \": \", substr(seq,i,i+w-1), sep=\"\" ); text( i-1, profile[i], h, pos=4, col=4 ); }\n"
			+ "}";

	/**
	 * This method creates a simple plot of the profile of scores for a sequence
	 * and a start position.
	 * 
	 * @param motifDisc
	 *            the {@link MotifDiscoverer}
	 * @param component
	 *            the index of the component
	 * @param motif
	 *            the index of the motif
	 * @param sequence
	 *            the given sequence
	 * @param startpos
	 *            the start position in the sequence
	 * @param r
	 *            the R environment that is used for drawing the plot
	 * @param width
	 *            the width of the image in pixel
	 * @param height
	 *            the height of the image in pixel
	 * @param kind
	 *            the kind of the score profile
	 * 
	 * @return an image packed in an {@link ImageResult}
	 * 
	 * @throws Exception
	 *             if something went wrong
	 */
	public static ImageResult plot(MotifDiscoverer motifDisc, int component,
			int motif, Sequence sequence, int startpos, REnvironment r,
			int width, int height, KindOfProfile kind) throws Exception {
		r.createVector("profile", motifDisc.getProfileOfScoresFor(component,
				motif, sequence, startpos, kind));
		r.voidEval("lim=c(min(profile),max(profile))");
		return new ImageResult("profile", "the profile of scores for "
				+ sequence.toString(startpos), r.plot(PLOT, width, height));
	}

	/**
	 * This method creates a plot of the profile of scores for a sequence and a
	 * start position and annotates bindings sites in the plot that have a higher
	 * score than <code>threshold</code>.
	 * 
	 * @param motifDisc
	 *            the {@link MotifDiscoverer}
	 * @param component
	 *            the index of the component
	 * @param motif
	 *            the index of the motif
	 * @param sequence
	 *            the given sequence
	 * @param startpos
	 *            the start position in the sequence
	 * @param r
	 *            the R environment that is used for drawing the plot
	 * @param width
	 *            the width of the image in pixel
	 * @param height
	 *            the height of the image in pixel
	 * @param yMin
	 *            the minimal value of the y-axis
	 * @param yMax
	 *            the maximal value of the y-axis
	 * @param threshold
	 *            the threshold for the annotation
	 * @param kind
	 *            the kind of the score profile
	 * 
	 * @return an image packed in a result
	 * 
	 * @throws Exception
	 *             if some thing went wrong
	 */
	public static ImageResult plotAndAnnotate(MotifDiscoverer motifDisc,
			int component, int motif, Sequence sequence, int startpos,
			REnvironment r, int width, int height, double yMin, double yMax,
			double threshold, KindOfProfile kind) throws Exception {
		r.createVector("profile", motifDisc.getProfileOfScoresFor(component,
				motif, sequence, startpos, kind));
		r.voidEval("lim=c(" + yMin + "," + yMax + ")");
		r.voidEval("threshold = " + threshold);
		String seq = sequence.toString(startpos);
		r.voidEval("seq = \"" + seq + "\"");
		r.voidEval("w = "
				+ motifDisc.getMotifLength(motifDisc
						.getGlobalIndexOfMotifInComponent(component, motif)));
		return new ImageResult("annotated profile with threshold " + threshold,
				"the annotated profile of scores for " + seq, r.plot(PLOT
						+ ANNOTATE, width, height));
	}
}
