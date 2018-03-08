package projects.gemoma.JunitTest;

import static org.junit.Assert.assertTrue;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Date;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import de.jstacs.parameters.ExpandableParameterSet;
import de.jstacs.parameters.FileParameter.FileRepresentation;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.ParameterSetContainer;
import de.jstacs.parameters.SelectionParameter;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.results.ResultSet;
import de.jstacs.results.TextResult;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ui.cli.CLI.SysProtocol;
import projects.gemoma.ExtractRNAseqEvidence;
import projects.gemoma.ExtractRNAseqEvidence.Stranded;
import projects.gemoma.Extractor;
import projects.gemoma.GeMoMa;
import projects.gemoma.GeMoMaAnnotationFilter;
import projects.gemoma.Tools.Ambiguity;

/**
 * 
 * @author Jens Keilwagen
 */
public class GeMoMaTest {

//for creating infrastructure, clean up, and nice things
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
		System.out.println("start: " + new Date());
	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
		System.out.println("end: " + new Date());
	}

	@Before
	public static void setUp() throws Exception {
		//this is a hack to avoid throwing useless exceptions while comparing files that include a version information
		/*
		Field f = GeMoMa.class.getField("version");
		f.setAccessible(true);
		
		Field modifiersField = Field.class.getDeclaredField("modifiers");
		modifiersField.setAccessible(true);
		modifiersField.setInt(f, f.getModifiers() & ~Modifier.FINAL);

		f.set(null, "0.0.0");//TODO
		/**/
	}

	@After
	public static void tearDown() throws Exception {
		System.out.println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
	}

//own test implementation	
	static String in = "projects/gemoma/JunitTest/given/";

	static void assertFile( String name, String givenName, ResultSet rs ) throws AssertionError, IOException {
		TextResult tr = (TextResult) (rs.getResultForName(name));
		FileRepresentation fr = tr.getValue();
		assertFile(name, in + givenName, fr.getFilename() );
	}
	
	static void assertFile( String name, String givenName, String fName ) throws AssertionError, IOException {	
		long l = -100;
		String l1, l2;
		File f1, f2;
		BufferedReader r1, r2;
		l1=l2=null;
		f1=f2=null;
		r1=r2=null;
		try {
			//File help = File.createTempFile(name.replaceAll(" ", "_"), ".tmp", new File(in));
			//FileManager.writeFile(help, fr.getContent() );
			f1 = new File(givenName);
			f2 = new File(fName);
			
			r1= new BufferedReader( new FileReader(f1) );
			r2= new BufferedReader( new FileReader(f2) );
			l=0;
			while( true ) {
				while( (l1=r1.readLine()) != null && l1.charAt(0)=='#' );//ignore comment lines
				while( (l2=r2.readLine()) != null && l2.charAt(0)=='#' );
				
				if( l1==null || l2 == null || !l1.equals(l2) ) {
					break;
				}
				l++;
			}
			r1.close();
			r2.close();
			if( l1 == null && l2 == null ) {
				//okay
			} else {
				throw new RuntimeException();
			}
		} catch( Exception e ) {
			if( l < 0 ) {
				throw new AssertionError((name==null?"":(name + ": "))+"The files differ."
						+ "\nFile 1 ("  + f1.getAbsolutePath() + "): " + (f1.exists() ? "exists" : "does not exist")
						+ "\nFile 2 ("  + f2.getAbsolutePath() + "): " + (f2.exists() ? "exists" : "does not exist"),
						e );
			} else {
				throw new AssertionError((name==null?"":(name + ": "))+"The files differ in line "+ l+"."
						+ "\nFile 1: "+f1.getAbsolutePath()+"\nContent: " + (l1==null?"[EOF]":l1) 
						+ "\nFile 2: "+f2.getAbsolutePath()+"\nContent: " + (l2==null?"[EOF]":l2) );
			}
		} finally {
			if( r1 != null ) {
				r1.close();
			}
			if( r2 != null ) {
				r2.close();
			}
		}
		
		//assertTrue( name + ": The files differ!", FileUtils.contentEquals( f1, f2 ));
	}
	
	static void assertStructure( ResultSet[] rs, int num, int... len ) {
		assertTrue( "Different number of ResultSets", rs.length == num );
		for( int i = 0; i < len.length; i++ ) {
			assertTrue( "Different number of Results in ResultSet " + i, rs[i].getNumberOfResults() == len[i] );
		}
	}
	
	void set( ParameterSet ps, int idx, int i, String value ) throws IllegalValueException, CloneNotSupportedException {
		ParameterSetContainer psc = (ParameterSetContainer) ps.getParameterAt(idx);
		ExpandableParameterSet eps = (ExpandableParameterSet) psc.getValue();
		while( i >= eps.getNumberOfParameters() ) {
			eps.addParameterToSet();
		}
		SimpleParameterSet sps = (SimpleParameterSet) ((ParameterSetContainer) eps.getParameterAt(i)).getValue();
		sps.getParameterAt(0).setValue(in + value);
	}
	
	@Test
	public void testExtractor() throws Exception {
		Extractor e = new Extractor(-1);
		
		//parameters
		ToolParameterSet ps = e.getToolParameters();
		ps.getParameterForName("annotation").setValue( in + "Arabidopsis_lyrata.v.1.0.31.chr.gff3");
		ps.getParameterForName("genome").setValue( in + "Arabidopsis_lyrata.v.1.0.31.dna.genome.fa");
		ps.getParameterForName("repair").setValue( true );
		ps.getParameterForName("full-length").setValue( false );
		ps.getParameterForName("Ambiguity").setValue( Ambiguity.AMBIGUOUS );
		ps.getParameterForName("proteins").setValue( true );
					
		//get results
		ResultSet[] rs = e.run(ps, new SysProtocol(), new ProgressUpdater(), 1).getValue();
		
		//compare that everything is correct
		assertStructure( rs, 1, 3 );
		assertFile("cds-parts", "cds-parts.fasta", rs[0] );
		assertFile("assignment", "assignment.tabular", rs[0] );
		assertFile("proteins", "proteins.fasta", rs[0] );
	}
	
	@Test
	public void testERE() throws Exception {
		ExtractRNAseqEvidence e = new ExtractRNAseqEvidence();
		
		//parameters
		ToolParameterSet ps = e.getToolParameters();
		ps.getParameterForName("Stranded").setValue( Stranded.FR_UNSTRANDED );
		set( ps, 1, 0, "RNAseq.bam");
		ps.getParameterForName("coverage output").setValue( true );

		//get results
		ResultSet[] rs = e.run(ps, new SysProtocol(), new ProgressUpdater(), 1).getValue();
		
		//compare that everything is correct
		assertStructure( rs, 1, 2 );
		assertFile("coverage", "coverage.bedgraph", rs[0] );
		assertFile("introns", "introns.gff", rs[0] );
	}
	
	@Test
	public void testGeMoMa() throws Exception {
		GeMoMa e = new GeMoMa(-1, 3600, 60*60*24*7);
		
		boolean simple = false;
		
		//parameters
		ToolParameterSet ps = e.getToolParameters();
		ps.getParameterForName("tblastn results").setValue( in + "tblastn.tabular");
		ps.getParameterForName("contig threshold").setValue( 0.4 );
		ps.getParameterForName("predictions").setValue(10);
		ps.getParameterForName("target genome").setValue( in + "TAIR10_chr_all.fas");
		ps.getParameterForName("assignment").setValue( in + "assignment.tabular");
		ps.getParameterForName("cds parts").setValue( in + "cds-parts.fasta");
		ps.getParameterForName("query proteins").setValue( in + "proteins.fasta");
		if( simple ) {
			set( ps, 4, 0, "at-introns.gff");
		} else {
			set( ps, 4, 0, "at-introns-1.gff");
			set( ps, 4, 1, "at-introns-2.gff");
		}
		ExpandableParameterSet eps = (ExpandableParameterSet) ps.getParameterAt(7).getValue();
		if( simple ) {
			SimpleParameterSet x = (SimpleParameterSet) ((ParameterSetContainer) eps.getParameterAt(0)).getValue();	
			SelectionParameter sp = (SelectionParameter) x.getParameterAt(0);
			sp.setValue("UNSTRANDED");
			((SimpleParameterSet) sp.getValue()).getParameterAt(0).setValue( in + "coverage.bedgraph");
		} else {
			SimpleParameterSet x = (SimpleParameterSet) ((ParameterSetContainer) eps.getParameterAt(0)).getValue();	
			SelectionParameter sp = (SelectionParameter) x.getParameterAt(0);
			sp.setValue("UNSTRANDED");
			((SimpleParameterSet) sp.getValue()).getParameterAt(0).setValue( in + "coverage.bedgraph-0.txt");
			
			eps.addParameterToSet();
			x = (SimpleParameterSet) ((ParameterSetContainer) eps.getParameterAt(1)).getValue();	
			sp = (SelectionParameter) x.getParameterAt(0);
			sp.setValue("UNSTRANDED");
			((SimpleParameterSet) sp.getValue()).getParameterAt(0).setValue( in + "coverage.bedgraph-1.txt");
		}
		
		//only first part
		//ps.getParameterForName("selected").setValue( in + "selected-first.txt" );
		
		//get results
		ResultSet[] rs = e.run(ps, new SysProtocol(), new ProgressUpdater(), 1).getValue();
		
		//compare that everything is correct
		assertStructure( rs, 1, 2 );
		assertFile("predicted annotation", "predicted_annotation.gff", rs[0] );
		assertFile("predicted protein", "predicted_protein.fasta", rs[0] );
	}
	
	@Test
	public void testGAF() throws Exception {
		GeMoMaAnnotationFilter e = new GeMoMaAnnotationFilter();
		
		//parameters
		ToolParameterSet ps = e.getToolParameters();
		set( ps, 6, 0, "predicted_annotation.gff");
		
		//get results
		ResultSet[] rs = e.run(ps, new SysProtocol(), new ProgressUpdater(), 1).getValue();
		
		//compare that everything is correct
		assertStructure( rs, 1, 1 );
		assertFile("filtered predictions", "filtered_predictions.gff", rs[0] );
	}
}
