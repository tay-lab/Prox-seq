import java.io.*;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.TreeSet;
import java.util.Objects;
import java.util.zip.*;

//import java.lang.Math;
import java.time.LocalDate;
import java.time.LocalTime;
import java.time.LocalDateTime;

import org.apache.commons.lang3.StringUtils;
import org.apache.commons.text.similarity.LevenshteinDistance;
import com.google.common.collect.*;

public class PLA_alignment_debug
{

	public static void main(String[] args)
	{
		switch (args[0])
		{
		
		/*
		 * Check if there are duplications across different cell barcodes for 2 modes: PLA product (PLA ID + UMI) vs UMI alone
		 * 		dirI=... MODE=... SUMMARY=...
		 * Input arguments:
		 * 		dirI: directory containing the .txt.gz files to process
		 * 		MODE: PLA or UMI, to check duplication of PLA products or UMI, respectively
		 * 		SUMMARY: path to saving the summary of UMI duplication
		 * NOTE: the files in dirI should have gone through UMI merging!!!
		 */
		case "Duplication":
		{
			// Parse the arguments
			String dirI = "", MODE = "", SUMMARY = "";
			for (int i=1; i<args.length; i++)
			{
				String[] j = args[i].split("=");
				
				switch (j[0])
				{
				case "dirI": dirI = j[1]; break;
				case "MODE": MODE = j[1]; break;
				case "SUMMARY": SUMMARY = j[1]; break;
				default: throw new java.lang.IllegalArgumentException("Invalid argument key specifier!");
				}
			}
			
			try (
					BufferedWriter bwsum = new BufferedWriter(new FileWriter(SUMMARY)) // summary file
				)
			{
				
				// Read the files in directory dirI
				File dir = new File(dirI);
				File[] dir_list = dir.listFiles();
				if (dir_list != null)
				{
					Multiset<String> UMIcount = HashMultiset.create(); // count the occurences of each UMI
					for (File child : dir_list)
					{
						try (
								BufferedReader brI = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(child))))
							)
						{
							
							// Read the PLA products and store the counts of each UMI in the HashMultiset
							String line;
							while ((line = brI.readLine()) != null)
							{
								
								switch (MODE)
								{
								case "PLA":
								{
									String[] values = line.split(",", 2);
									UMIcount.add(values[1]);
									break;
								}
								case "UMI":
								{
									String[] values = line.split(",");
									UMIcount.add(values[1]);
									break;
								}
								default: throw new java.lang.IllegalArgumentException("Invalid MODE argument!");
								}
								
								
							}
							
						} catch (IOException e) {throw new IllegalArgumentException("Invalid file format in the input directory!");}
						
					}
					
					// Sort the HashMultiset by decreasing occurences, and save to an array
					String[] UMIcount_sorted = Multisets.copyHighestCountFirst(UMIcount).elementSet().toArray(new String[0]);
					
					// Print out the UMI and counts in SUMMARY file
					bwsum.write("UMI\tcount");
					bwsum.newLine();
					for (int j=0; j<UMIcount_sorted.length; j++)
					{
						bwsum.write(UMIcount_sorted[j] + "\t" + UMIcount.count(UMIcount_sorted[j]));
						bwsum.newLine();
					}
				}
				else
				{
					throw new IllegalArgumentException("Empty file directory!");
				}
				
			} catch (IOException e) {throw new IllegalArgumentException("Invalid file paths!");}
			
			break;
		}
		
		/*
		 * Save reads that have a specific AB barcode
		 * 		R1=... O=... AB=...
		 * Input arguments:
		 * 		R1: path to read 1 file (fastq.gz format)
		 * 		O: path to save raw reads (txt.gz format)
		 * 		AB: the sequence of the AB barcode to look for (eg, AGAGTCTC)
		 */
		
		case "MatchAB":
		{
			// Parse arguments
			String R1 = "", O = "", AB = "";
			for (int i=1; i<args.length; i++)
			{
				String[] j = args[i].split("=");
				
				switch (j[0])
				{
				case "R1": R1 = j[1]; break;
				case "O": O = j[1]; break;
				case "AB": AB = j[1]; break;
				default: throw new java.lang.IllegalArgumentException("Invalid argument key specifier!");
				}
			}
			
			try (
					BufferedReader br1 = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(R1)))); // read1
					BufferedWriter bwout = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(O)))) // output file
				)
			{
				
				// Declare variables
				String line1;
				int counter = 0;
				String connector = "TCGTGTCGTGTCGTGTCTAAAG"; // connector sequence	
				
				bwout.write("Raw read\tFrame shift\tConnector region\tAB1\tAB2"); bwout.newLine();
				
				while ((line1=br1.readLine()) != null)
				{
							
					if ((counter % 4) == 1)
					{				
						
						// Check if the UMI region has an excessive amount of A
						if (StringUtils.countMatches(line1.substring(1, 17), "A") >= 10)
						{
							counter++;
							continue;
						}
						
						// Expected starting location of connector at index 39 (base 40th)
						int connector_start = 39;
						int connector_start_temp = 39; // temporary starting location, used in for loop
						
						// Locate the connector using Levenshtein distance, max allowed distance is 2
						// Does not need to take into account N, since the connector region doesn't usually contain N
						int[] connector_shift = new int[] {0, -1, 1, -2, 2, -3, 3};
						boolean match_connector = false; // true if matching connector is found
						int connector_distance = 100; // lowest levenshtein distance found between the true connector sequence and the test connector sequence
						for (int shift_i : connector_shift)
						{
							int temp_distance = LevenshteinDistance.getDefaultInstance().apply(line1.substring(connector_start+shift_i, connector_start+shift_i+connector.length()), connector);
							if (temp_distance <= 2)
							{
								if (temp_distance < connector_distance) // update the connector_distance if found a better matching test connector sequence
								{
									connector_distance = temp_distance;
									connector_start_temp = connector_start + shift_i;
									match_connector = true;
								}
								if (temp_distance == 0)
								{break;}
							}
						}
						connector_start = connector_start_temp;
						
						if (match_connector)
						{
							// Check if there is a frameshift
							int shift_j = line1.substring(connector_start, connector_start+connector.length()+1).indexOf("TAAAG"); // location of TAAAG in the found connector region
							if (shift_j == -1)
							{
								counter++;
								continue;
							}
							shift_j = 17 - shift_j; // number of bases to shift to the left
							
							// Found AB barcodes
							String AB1_found = line1.substring(connector_start-18, connector_start-18+8);
							String AB2_found = line1.substring(connector_start+25-shift_j, connector_start+25-shift_j+8);
			
							if ((Objects.equals(AB1_found, AB)) || (Objects.equals(AB2_found, AB)))
							{
								bwout.write(line1 + "\t" + String.format("%d", shift_j) + "\t" + line1.substring(connector_start, connector_start+connector.length()) +
										"\t" + AB1_found + "\t" + AB2_found); bwout.newLine();
							}
						}
						
						
						
					}
					
					
//					if ((((counter-1)/4)+1)>2000000) {System.out.printf("%s   %s   Processed %,d lines%n", LocalDate.now(), LocalTime.now().withNano(0), counter); break;} // for testing purposes
				
					counter++;
				}
				
			} catch (IOException e) {throw new IllegalArgumentException("Invalid file paths!");}
			
			
			break;
		}
		
		/*
		 * Save reads that returns OutOfIndex error
		 * 		R1=... O=...
		 * Input arguments:
		 * 		R1: path to read 1 file (fastq.gz format)
		 * 		O: path to save raw reads (txt.gz format)
		 * 		AB: the sequence of the AB barcode to look for (eg, AGAGTCTC)
		 */
		
		case "OutOfBound":
		{
			// Parse arguments
			String R1 = "", O = "";
			for (int i=1; i<args.length; i++)
			{
				String[] j = args[i].split("=");
				
				switch (j[0])
				{
				case "R1": R1 = j[1]; break;
				case "O": O = j[1]; break;
				default: throw new java.lang.IllegalArgumentException("Invalid argument key specifier!");
				}
			}
			
			try (
					BufferedReader br1 = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(R1)))); // read1
					BufferedWriter bwout = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(O)))) // output file
				)
			{
				// Declare variables
				String line1;
				int counter = 0;
				String connector = "TCGTGTCGTGTCGTGTCTAAAG"; // connector sequence	
				
				while ((line1=br1.readLine()) != null)
				{
							
					if ((counter % 4) == 1)
					{				
//						bwout.write(line1);
//						bwout.newLine();
						
						// Check if the UMI region has an excessive amount of A
						if (StringUtils.countMatches(line1.substring(1, 17), "A") >= 10)
						{
							counter++;
							continue;
						}
						
						// Expected starting location of connector at index 39 (base 40th)
						int connector_start = 39;
						int connector_start_temp = 39; // temporary starting location, used in for loop
						
						// Locate the connector using Levenshtein distance, max allowed distance is 2
						// Does not need to take into account N, since the connector region doesn't usually contain N
						int[] connector_shift = new int[] {0, -1, 1, -2, 2, -3, 3, -4, 4};
						boolean match_connector = false; // true if matching connector is found
						int connector_distance = 100; // lowest levenshtein distance found between the true connector sequence and the test connector sequence
						for (int shift_i : connector_shift)
						{
							int temp_distance = LevenshteinDistance.getDefaultInstance().apply(line1.substring(connector_start+shift_i, connector_start+shift_i+connector.length()), connector);
							if (temp_distance <= 2)
							{
								if (temp_distance < connector_distance) // update the connector_distance if found a better matching test connector sequence
								{
									connector_distance = temp_distance;
									connector_start_temp = connector_start + shift_i;
									match_connector = true;
								}
								if (temp_distance == 0)
								{break;}
							}
						}
						connector_start = connector_start_temp;
						
						if (match_connector)
						{
							// Check if there is a frameshift
							int shift_j = line1.substring(connector_start, connector_start+connector.length()+1).indexOf("TAAAG"); // location of TAAAG in the found connector region
							if (shift_j == -1)
							{
								counter++;
								continue;
							}
							shift_j = 17 - shift_j; // number of bases to shift to the left
							
							// Found AB barcodes
							if ((connector_start+25-shift_j+8) > 75)
							{
								bwout.write(line1);
								bwout.newLine();
							}

			

						}
					}
					counter++;
				}
					
			} catch (IOException e) {throw new IllegalArgumentException("Invalid file paths!");}
			
			break;
		}
		
		default: throw new java.lang.IllegalArgumentException("Invalid function call!");
			
		}
		
	
		

	}

}
