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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URISyntaxException;
import java.net.URL;
import java.nio.channels.Channels;
import java.nio.channels.ReadableByteChannel;
import java.util.ArrayList;

import de.jstacs.DataType;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.ParameterException;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;

public class NCBIReferenceRetriever extends GeMoMaModule {

	@Override
	public ToolParameterSet getToolParameters() {
		try {
			return new ToolParameterSet( getToolName(), 
					new SimpleParameter(DataType.STRING, "reference directory", "the directory where the genome and annotation files of the reference organisms should be stored", true, "references/"),
					new SimpleParameter(DataType.INT, "number of tries", "the number of tries for downloading a reference file", true, new NumberValidator<Integer>(1, 100), 10),
					new FileParameter("reference list", "a list of reference organisms", "txt", true)
			);
		} catch (ParameterException e) {
			e.printStackTrace();
			throw new RuntimeException();
		}
	}

	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads, String temp)
			throws Exception {
		String base = "https://www.ncbi.nlm.nih.gov/genome/?term=";
		
		//user-parameters
		String outdir = (String) parameters.getParameterForName("reference directory").getValue();
		int max = (Integer) parameters.getParameterForName("number of tries").getValue();
		FileParameter fp = (FileParameter) parameters.getParameterForName("reference list");
		
		//read reference file
		ArrayList<String> ref = new ArrayList<String>();
		
		BufferedReader re = new BufferedReader( new FileReader( fp.getValue().toString() ) );
		String line;
		while( (line=re.readLine()) != null ) {
			ref.add(line);
		}
		re.close();
		//{ "arabidopsis thaliana", "mickey mouse"};
		
		
		progress.setLast(ref.size()*2);
		progress.setCurrent(0);
		
		BufferedReader buffer = null;
		ArrayList<String> needToBeDownloaded = new ArrayList<String>(), 
				downloaded = new ArrayList<String>(),
				check = new ArrayList<String>();
		
		for( String r: ref ) {
			String url = base+r.replaceAll(" ", "+");
			protocol.append( r + "\t" + url + "\n" );
			
			buffer = new BufferedReader(new InputStreamReader(new URL(url).openStream(), "iso-8859-9"));
			StringBuffer builder = new StringBuffer();
			while ((line=buffer.readLine()) != null ) {
				builder.append( line + "\n");
			}
			buffer.close();
			
			String genome = getLink(builder, "Download sequences in FASTA format for");
			String annotation = getLink(builder, "Download genome annotation in");
			
			if( genome.length()==0 || annotation.length()==0 ) {
				check.add( r );
			}
			
			//System.out.println( "\t" + genome + "\t" + annotation );
			
			download( protocol, max, outdir, genome, downloaded, needToBeDownloaded );
			progress.add(1);
			download( protocol, max, outdir, annotation, downloaded, needToBeDownloaded );
			progress.add(1);
			protocol.append( "\n" );
		}
		
		if( downloaded.size()>0 ) {
			protocol.append( "Downloaded " + downloaded.size() + " file(s):\n" );
			for( String f : downloaded ) {
				protocol.append( f + "\n" );
			}
		}

		protocol.append( "\n" );
		if( needToBeDownloaded.size()==0 ) {
			protocol.append("Nothing to download\n");
		} else {
			protocol.append("Please download " + needToBeDownloaded.size() + " file(s):\n");
			for( String f : needToBeDownloaded ) {
				protocol.append(f);
			}
		}
		
		if( check.size()>0 ) {
			protocol.append( "\n" );
			protocol.appendWarning("Please check the following " + check.size() + " reference species:\n");
			for( String f : check ) {
				protocol.appendWarning( f + "\n" );
			}
		}
		
		//TODO
		return null;
	}
	
	static String getLink( StringBuffer content, String pattern ) {
		int idx = content.indexOf(pattern);
		if( idx>=0 ) {
			idx+=pattern.length();
			idx=content.indexOf("href=\"",idx)+6;
			return content.substring(idx,content.indexOf("\"",idx));
		} else {
			return "";
		}
	}
	
	static void download( Protocol protocol, int max, String outdir, String link, ArrayList<String> downloaded, ArrayList<String> needToBeDownloaded ) throws URISyntaxException, IOException {
		if( link.length()!=0 ) {
			int idx = link.lastIndexOf("/");
			String check = outdir + link.substring(idx+1);
			
			File f = new File( check );
			if( !f.exists() ) {
protocol.append(link+"\n");
				boolean successful = false;
				Exception myE = null;
				URL website = new URL(link);
				int i = 1;
				do {
protocol.append("try " + i++ + "\n");
					try {
						ReadableByteChannel rbc = Channels.newChannel(website.openStream());
						FileOutputStream fos = new FileOutputStream(outdir + link.substring(link.lastIndexOf("/")) );
						fos.getChannel().transferFrom(rbc, 0, Long.MAX_VALUE);
						fos.close();
						rbc.close();
						successful = true;
					} catch( Exception e ) {
						myE = e;
					}
					
					if( !successful ) {
						if( !myE.getMessage().startsWith("sun.net.ftp.FtpProtocolException") ) {
							protocol.appendThrowable(myE);
							break;
						}
					} else {
						break;
					}
					
				} while( i<max );
		        if( successful ) {
		        	downloaded.add(link);
		        } else {
					needToBeDownloaded.add(link);
		        }
			}
		}
	}

	@Override
	public String getToolName() {
		return "NCBI Reference Retriever";
	}

	@Override
	public String getShortName() {
		return "NRR";
	}

	@Override
	public String getDescription() {
		return "downloads new assembly and annotation files of reference organisms";
	}

	@Override
	public String getHelpText() {
		return "This tool can be used to download or update assembly and annotation files of reference organsims from NCBI."
				+ " This way it allows to easily collect all data necessary to start **GeMoMaPipeline** or **Extractor**."
				+ MORE;
	}

	@Override
	public ResultEntry[] getDefaultResultInfos() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public ToolResult[] getTestCases( String path ) {
		// TODO Auto-generated method stub
		return null;
	}
}
