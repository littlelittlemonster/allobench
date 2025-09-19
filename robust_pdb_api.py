#!/usr/bin/env python3
"""
Robust PDB API functions with comprehensive error handling
"""

import requests
import json
import time
import pandas as pd
from typing import List, Dict, Any
import random

class RobustPDBAPI:
    """Robust PDB API client with retry logic and error handling"""
    
    def __init__(self, max_retries=5, base_delay=1, max_delay=60):
        self.max_retries = max_retries
        self.base_delay = base_delay
        self.max_delay = max_delay
        self.session = requests.Session()
        
        # Set headers to mimic a real browser
        self.session.headers.update({
            'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36',
            'Accept': 'application/json',
            'Accept-Language': 'en-US,en;q=0.9',
            'Accept-Encoding': 'gzip, deflate, br',
            'Connection': 'keep-alive',
            'Upgrade-Insecure-Requests': '1',
        })
    
    def exponential_backoff(self, attempt: int) -> float:
        """Calculate exponential backoff delay with jitter"""
        delay = min(self.base_delay * (2 ** attempt), self.max_delay)
        jitter = random.uniform(0.1, 0.3) * delay
        return delay + jitter
    
    def make_request(self, query: str, variables: Dict[str, Any], timeout: int = 120) -> Dict[str, Any]:
        """Make a robust request to PDB GraphQL API"""
        
        for attempt in range(self.max_retries):
            try:
                response = self.session.post(
                    url='https://data.rcsb.org/graphql',
                    json={'query': query, 'variables': variables},
                    timeout=timeout
                )
                
                if response.status_code == 200:
                    data = response.json()
                    if 'errors' in data:
                        print(f"GraphQL errors: {data['errors']}")
                        return None
                    return data
                elif response.status_code == 429:  # Rate limited
                    print(f"Rate limited (429), waiting longer...")
                    wait_time = self.exponential_backoff(attempt) * 2
                    time.sleep(wait_time)
                    continue
                else:
                    print(f"HTTP error {response.status_code}: {response.text[:200]}")
                    
            except requests.exceptions.ConnectionError as e:
                print(f"Connection error (attempt {attempt + 1}/{self.max_retries}): {e}")
            except requests.exceptions.Timeout as e:
                print(f"Timeout error (attempt {attempt + 1}/{self.max_retries}): {e}")
            except requests.exceptions.RequestException as e:
                print(f"Request error (attempt {attempt + 1}/{self.max_retries}): {e}")
            except Exception as e:
                print(f"Unexpected error (attempt {attempt + 1}/{self.max_retries}): {e}")
            
            if attempt < self.max_retries - 1:
                wait_time = self.exponential_backoff(attempt)
                print(f"Waiting {wait_time:.1f} seconds before retry...")
                time.sleep(wait_time)
        
        print(f"Failed after {self.max_retries} attempts")
        return None
    
    def get_pdb_data(self, pdb_ids: List[str], batch_size: int = 20) -> pd.DataFrame:
        """Get PDB data with robust error handling"""
        
        query = """query structure($pdb_ids: [String!]!) {
  entries(entry_ids: $pdb_ids) {
    rcsb_id
    rcsb_entry_info {
          experimental_method
          resolution_combined
        }
    assemblies {
      rcsb_id
      polymer_entity_instances {
        rcsb_id
        polymer_entity {
          entity_poly {
            pdbx_seq_one_letter_code_can
          }
          uniprots {
            rcsb_id
          }
        }
      }
      rcsb_struct_symmetry {
        kind
        oligomeric_state
        stoichiometry
      }
    }
  }
}"""
        
        all_pdb_data = []
        total_ids = len(pdb_ids)
        successful_batches = 0
        failed_batches = 0
        
        print(f"Processing {total_ids} PDB IDs in batches of {batch_size}")
        
        for i in range(0, total_ids, batch_size):
            batch_ids = pdb_ids[i:i + batch_size]
            batch_num = i // batch_size + 1
            total_batches = (total_ids + batch_size - 1) // batch_size
            
            print(f"Processing batch {batch_num}/{total_batches} ({len(batch_ids)} IDs)")
            
            data = self.make_request(query, {'pdb_ids': list(batch_ids)})
            
            if data and 'data' in data and data['data']:
                entries = data['data']['entries']
                print(f"  ✓ Successfully retrieved {len(entries)} entries")
                
                for entry in entries:
                    pdb_id = entry['rcsb_id']
                    experimental_method = entry['rcsb_entry_info']['experimental_method']
                    resolution = entry['rcsb_entry_info']['resolution_combined']
                    
                    if resolution and len(resolution) == 1:
                        resolution = resolution[0]
                    
                    first_assembly = entry['assemblies'][0] if entry['assemblies'] else None
                    
                    chain_uniprot_mapping = []
                    if first_assembly:
                        for polymer_entity in first_assembly['polymer_entity_instances']:
                            chain = polymer_entity['rcsb_id'].split('.')[1]
                            sequence = polymer_entity['polymer_entity']['entity_poly']['pdbx_seq_one_letter_code_can']
                            polymer_entity_uniprot = polymer_entity['polymer_entity']['uniprots']
                            uniprot_ids = []
                            if polymer_entity_uniprot:
                                for uniprot_id in polymer_entity_uniprot:
                                    uniprot_ids.append(uniprot_id['rcsb_id'])
                            chain_uniprot_mapping.append([chain, uniprot_ids, sequence])
                    
                    oligomeric_state = None
                    stoichiometry = None
                    if first_assembly and first_assembly['rcsb_struct_symmetry'] is not None:
                        for symmetry in first_assembly['rcsb_struct_symmetry']:
                            if symmetry['kind'] == 'Global Symmetry':
                                oligomeric_state = symmetry['oligomeric_state']
                                stoichiometry = symmetry['stoichiometry']
                    
                    all_pdb_data.append([
                        pdb_id, chain_uniprot_mapping, oligomeric_state, 
                        stoichiometry, experimental_method, resolution
                    ])
                
                successful_batches += 1
            else:
                print(f"  ✗ Failed to retrieve data for batch {batch_num}")
                failed_batches += 1
            
            # Add delay between batches
            if i + batch_size < total_ids:
                time.sleep(2)
        
        print(f"\nSummary:")
        print(f"  Successful batches: {successful_batches}/{total_batches}")
        print(f"  Failed batches: {failed_batches}/{total_batches}")
        print(f"  Total entries retrieved: {len(all_pdb_data)}")
        
        return pd.DataFrame(all_pdb_data, columns=[
            'PDB ID', 'Map PDB Chain to UniProt', 'Oligomeric State', 
            'Stoichiometry', 'Experimental Method', 'Resolution'
        ])

# Convenience function for backward compatibility
def get_pdb_data_robust(pdb_ids, batch_size=20):
    """Robust version of get_pdb_data function"""
    api = RobustPDBAPI()
    return api.get_pdb_data(pdb_ids, batch_size)

if __name__ == "__main__":
    # Test the robust function
    import sys
    import os
    sys.path.append('src')
    import parse_asd_xml
    
    # Load a small sample
    xml_dir = '../xml_file'
    df_asd = parse_asd_xml.asd_to_df(xml_dir)
    test_pdb_ids = df_asd['PDB ID'].unique()[:10]
    
    print(f"Testing robust PDB API with {len(test_pdb_ids)} PDB IDs")
    result = get_pdb_data_robust(test_pdb_ids, batch_size=5)
    print(f"Successfully retrieved data for {len(result)} PDB IDs")
