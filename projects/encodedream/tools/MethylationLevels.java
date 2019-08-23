package projects.encodedream.tools;
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
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Date;
import java.util.LinkedList;
import java.util.zip.GZIPOutputStream;

import de.jstacs.DataType;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.SelectionParameter;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameter.DatatypeNotValidException;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.results.ResultSet;
import de.jstacs.results.TextResult;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;
import projects.encodedream.Coverage;
import projects.encodedream.ObjectStream;
import projects.encodedream.Pileup;
import projects.encodedream.Pileup.CovPile;

public class MethylationLevels implements JstacsTool {

	@Override
	public ToolParameterSet getToolParameters() {
		LinkedList<Parameter> pars = new LinkedList<>();
		
		try {
			
			pars.add( new FileParameter("Input Bed.gz", "The bedMethyl file (gzipped) containing the methylation levels", "bed.gz", true) );
			pars.add( new FileParameter("FastA index", "The genome index", "fai", true) );
		} catch (Exception e1) {
			e1.printStackTrace();
		}

		try {
			pars.add(new SimpleParameter(DataType.INT, "Bin width", "The width of the genomic bins considered", true));
		} catch (DatatypeNotValidException e) {
			e.printStackTrace();
		}
		return new ToolParameterSet(getShortName(),pars.toArray(new Parameter[0]));
	}

	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads)
			throws Exception {
		
		File out = File.createTempFile("methyl", ".temp.gz", new File("."));
		out.deleteOnExit();

		GZIPOutputStream os = new GZIPOutputStream(new FileOutputStream(out));


		String methFile = ((FileParameter)parameters.getParameterAt(0)).getFileContents().getFilename();
		int bin = (int) parameters.getParameterAt(2).getValue();

		String faiFile = ((FileParameter)parameters.getParameterAt(1)).getFileContents().getFilename();
			
		Coverage.coverageMeth(faiFile, methFile, new PrintStream(os), bin);
		
		os.close();
		
		TextResult tr = new TextResult("Methylation levels", "Features computed from methylation levels", new FileParameter.FileRepresentation(out.getAbsolutePath()), "tsv.gz", getToolName(), null, true);
		
		return new ToolResult("Result of "+getToolName(), getToolName(), null, new ResultSet(tr), parameters, getToolName(), new Date(System.currentTimeMillis()) );
		
	}

	@Override
	public String getToolName() {
		return "Methylation levels";
	}

	@Override
	public String getToolVersion() {
		return "0.1";
	}

	@Override
	public String getShortName() {
		return "methyl";
	}

	@Override
	public String getDescription() {
		return "computes methylation level features";
	}

	@Override
	public String getHelpText() {
		return "";//TODO FIXME
	}

	@Override
	public ResultEntry[] getDefaultResultInfos() {
		return null;
	}


	@Override
	public ToolResult[] getTestCases(String path) {
		return null;
	}

	@Override
	public void clear() {		
	}

	@Override
	public String[] getReferences() {
		return null;
	}
	
}
