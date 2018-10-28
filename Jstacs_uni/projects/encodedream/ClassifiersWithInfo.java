package projects.encodedream;
/*
 * This file is part of Jstacs.
 * 
 * Jstacs is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 * 
 * Jstacs is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * Jstacs. If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */
import de.jstacs.Storable;
import de.jstacs.classifiers.differentiableSequenceScoreBased.gendismix.GenDisMixClassifier;
import de.jstacs.io.NonParsableException;
import de.jstacs.io.XMLParser;

public class ClassifiersWithInfo implements Storable {

	private GenDisMixClassifier[] classifiers;
	private int numBins;
	private int binWidth;
	private int numMotifs;
	private int binsBefore;
	private int binsAfter;
	
	public ClassifiersWithInfo(GenDisMixClassifier[] classifiers, int numBins, int binWidth, int binsBefore, int binsAfter, int numMotifs){
		this.classifiers = classifiers;
		this.numBins = numBins;
		this.binWidth = binWidth;
		this.numMotifs = numMotifs;
		this.binsBefore = binsBefore;
		this.binsAfter = binsAfter;
	}
	
	public ClassifiersWithInfo(StringBuffer xml) throws NonParsableException{
		classifiers = (GenDisMixClassifier[]) XMLParser.extractObjectForTags(xml, "classifiers");
		numBins = (int) XMLParser.extractObjectForTags(xml, "numBins");
		binWidth = (int) XMLParser.extractObjectForTags(xml, "binWidth");
		numMotifs = (int) XMLParser.extractObjectForTags(xml, "numMotifs");
		binsBefore = (int) XMLParser.extractObjectForTags(xml, "binsBefore");
		binsAfter = (int) XMLParser.extractObjectForTags(xml, "binsAfter");
	}
	
	public StringBuffer toXML(){
		StringBuffer xml = new StringBuffer();
		XMLParser.appendObjectWithTags(xml, classifiers, "classifiers");
		XMLParser.appendObjectWithTags(xml, numBins, "numBins");
		XMLParser.appendObjectWithTags(xml, binWidth, "binWidth");
		XMLParser.appendObjectWithTags(xml, numMotifs, "numMotifs");
		XMLParser.appendObjectWithTags(xml, binsBefore, "binsBefore");
		XMLParser.appendObjectWithTags(xml, binsAfter, "binsAfter");
		return xml;
	}

	public int getBinsBefore() {
		return binsBefore;
	}

	public int getBinsAfter() {
		return binsAfter;
	}

	public GenDisMixClassifier[] getClassifiers() {
		return classifiers;
	}

	public int getNumBins() {
		return numBins;
	}

	public int getBinWidth() {
		return binWidth;
	}

	public int getNumberOfMotifs() {
		return numMotifs;
	}

	public void limitClassifiers(Integer numClass) {
		if(numClass > 0){
			GenDisMixClassifier[] temp = new GenDisMixClassifier[numClass];
			for(int i=0;i<temp.length;i++){
				temp[i] = classifiers[i];
			}
			this.classifiers = temp;
		}else{
			GenDisMixClassifier[] temp = new GenDisMixClassifier[-numClass];
			for(int i=0;i<temp.length;i++){
				temp[i] = classifiers[classifiers.length-i-1];
			}
			this.classifiers = temp;
		}
	}
	
}
