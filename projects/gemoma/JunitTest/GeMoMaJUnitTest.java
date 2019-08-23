package projects.gemoma.JunitTest;

import static org.junit.Assert.assertTrue;

import java.util.ArrayList;

import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolResult;
import de.jstacs.tools.ui.cli.CLI.SysProtocol;
import projects.gemoma.AnnotationEvidence;
import projects.gemoma.AnnotationFinalizer;
import projects.gemoma.CheckIntrons;
import projects.gemoma.CompareTranscripts;
import projects.gemoma.Denoise;
import projects.gemoma.ExtractRNAseqEvidence;
import projects.gemoma.Extractor;
import projects.gemoma.GeMoMa;
import projects.gemoma.GeMoMaAnnotationFilter;
import projects.gemoma.GeMoMaPipeline;

/**
 * Parameterized JUnitTest for GeMoMa modules.
 * 
 * @author Jens Keilwagen
 */
@RunWith(Parameterized.class)
public class GeMoMaJUnitTest {
	
	private JstacsTool tool;
	private ToolResult tr;
	
	public GeMoMaJUnitTest(JstacsTool tool, ToolResult tr) {
		this.tool = tool;
		this.tr = tr;
	}
	
	@Parameters(name="{index}: {0}")//TODO improve name
	public static ArrayList<Object[]> testCases() {
		int maxSize = -1;
		long timeOut=3600, maxTimeOut=60*60*24*7;
		
		JstacsTool[] jt = {//TODO extend?
				new GeMoMaPipeline(),
				
				new ExtractRNAseqEvidence(),
				new CheckIntrons(),
				new Denoise(),
				
				new Extractor(maxSize),
				new GeMoMa(maxSize, timeOut, maxTimeOut),
				new GeMoMaAnnotationFilter(),
				new AnnotationFinalizer(),
				new AnnotationEvidence(),
				new CompareTranscripts(),
		};
		
		ArrayList<Object[]> tests = new ArrayList<Object[]>();
		for( JstacsTool t: jt ) {
			ToolResult[] tr = t.getTestCases(".");
			if( tr!= null ) {
				for( int i = 0; i < tr.length; i++ ) {
					tests.add( new Object[] {t,tr[i]} );
				}
			}
		}
		return tests;
	}
	
	private static Protocol protocol = new SysProtocol();
	private static ProgressUpdater progress = new ProgressUpdater();
	
 	@Test
	public void test() {
 		for( int i = 0; i < 150; i++ ) {
 			protocol.append( "=" );
 		}
 		protocol.append( "\n" );
 		protocol.append( "Test for " + tr.getToolName() + "\n" );
 		
 		boolean eq = false;
 		try {
 			tool.clear();
 			eq = tr.equals( tool.run( tr.getToolParameters(), protocol, progress, 1 ) );
 		} catch( Exception e ) {
 			e.printStackTrace();
 		}
 		assertTrue( eq );
	}
}
