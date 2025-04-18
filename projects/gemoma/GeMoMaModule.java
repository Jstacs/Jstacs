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

package projects.gemoma;

import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;

/**
 * This is the abstract class of all {@link GeMoMa} modules that are available as {@link JstacsTool}. 
 * 
 * @author Jens Keilwagen
 * 
 * @see JstacsTool
 * @see GeMoMa
 */
public abstract class GeMoMaModule implements JstacsTool {

	public static final String VERSION = "1.9";
	
	public static final String INFO = "#SOFTWARE INFO: ";

	
	public static String MORE = "\n\nFor more information please visit http://www.jstacs.de/index.php/GeMoMa\n"
			+ "If you have any questions, comments or bugs, please check FAQs on our homepage, our github page https://github.com/Jstacs/Jstacs/labels/GeMoMa or contact jens.keilwagen@julius-kuehn.de";

	public static String[] REF = {
			"@article{Keilwagen:2016:GeMoMa,\n"
		    +" author = {Keilwagen, Jens and Wenk, Michael and Erickson, Jessica L. and Schattat, Martin H. and Grau, Jan and Hartung, Frank},\n"
		    +" title = {{Using intron position conservation for homology-based gene prediction}},\n"
		    +" journal = {Nucleic Acids Research},\n"
		    +" volume = {44},\n"
		    +" number = {9},\n"
		    +" pages = {e89-e89},\n"
		    +" year = {2016},\n"
		    +" month = {02},\n"
		    +" issn = {0305-1048},\n"
		    +" doi = {10.1093/nar/gkw092}\n"
		+"}\n",
		"@article{Keilwagen:2018:GeMoMa_RNAseq,\n"
			+" author = {Keilwagen, Jens and Hartung, Frank and Paulini, Michael and Twardziok, Sven O. and Grau, Jan},\n"
			+" title = {Combining RNA-seq data and homology-based gene prediction for plants, animals and fungi},\n"
			+" journal = {BMC Bioinformatics},\n"
			+" year = {2018},\n"
			+" month = {May},\n"
			+" day = {30},\n"
			+" volume = {19},\n"
			+" number = {1},\n"
			+" pages = {189},\n"
			+" issn = {1471-2105},\n"
			+" doi = {10.1186/s12859-018-2203-5}\n"
		+"}\n"
	};
		
	@Override
	public final String[] getReferences() {
		return REF;
	}
	
	@Override
	public final String getToolVersion() {
		return VERSION;
	}
	
	public synchronized void clear() {
		GeMoMa.seqs=null;
		GeMoMa.selected=null;
		GeMoMa.donorSites=null;
		GeMoMa.acceptorSites=null;
		GeMoMa.coverage=null;
	}
	
	public final ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads) throws Exception {
		return run(parameters, protocol, progress, threads, Tools.GeMoMa_TEMP);
	}

	public abstract ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads, String temp ) throws Exception;
	
	/**
	 * Returns the requested parameter from user-given {@link ToolParameterSet}.
	 * If the user-given {@link ToolParameterSet} does not contain such a parameter, it is retrieved from the default {@link ToolParameterSet}.
	 * This method is especially useful for test cases if a new {@link Parameter} is introduced to a new version of a tool.
	 * 
	 * @param given the user given parameter set
	 * @param parameterName the name of the requested parameter
	 * 
	 * @return the requested parameter
	 * 
	 * @throws Exception if the default parameters cannot be created
	 * 
	 * @see #getToolParameters()
	 * @see #getTestCases(String)
	 */
	public Parameter getParameter( ToolParameterSet given, String parameterName ) throws Exception {
		Parameter p = given.getParameterForName( parameterName );
		if( p == null ) {
			p = getToolParameters().getParameterForName( parameterName );
		}
		return p;
	}
}
