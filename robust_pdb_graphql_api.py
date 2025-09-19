#!/usr/bin/env python3
"""
Robust PDB GraphQL API with error handling and retry logic
"""

import requests
import pandas as pd
import time
from tqdm import tqdm
import sys
from typing import List, Set
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

class RobustPDBGraphQLAPI:
    """Robust PDB GraphQL API with comprehensive error handling"""
    
    def __init__(self, max_retries=5, base_delay=2, max_delay=60):
        self.max_retries = max_retries
        self.base_delay = base_delay
        self.max_delay = max_delay
        
        # Create a session with retry strategy
        self.session = requests.Session()
        
        # Configure retry strategy
        retry_strategy = Retry(
            total=max_retries,
            backoff_factor=base_delay,
            status_forcelist=[429, 500, 502, 503, 504],
            allowed_methods=["HEAD", "GET", "POST", "OPTIONS"]
        )
        
        adapter = HTTPAdapter(max_retries=retry_strategy)
        self.session.mount("http://", adapter)
        self.session.mount("https://", adapter)
        
        # Set headers
        self.session.headers.update({
            'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36',
            'Accept': 'application/json',
            'Content-Type': 'application/json'
        })
    
    def exponential_backoff(self, attempt: int) -> float:
        """Calculate exponential backoff delay"""
        delay = min(self.base_delay * (2 ** attempt), self.max_delay)
        return delay
    
    def make_graphql_request(self, query: str, timeout: int = 120):
        """Make a robust GraphQL request with retry logic"""
        
        url = 'https://data.rcsb.org/graphql'
        
        for attempt in range(self.max_retries):
            try:
                print(f"  Making GraphQL request (attempt {attempt + 1}/{self.max_retries})...")
                
                response = self.session.post(
                    url=url,
                    json={"query": query},
                    timeout=timeout
                )
                
                # Check for HTTP errors
                response.raise_for_status()
                
                # Parse JSON response
                data = response.json()
                
                # Check for GraphQL errors
                if 'errors' in data:
                    print(f"  → GraphQL errors: {data['errors']}")
                    return None
                
                print(f"  ✓ Successfully retrieved data")
                return data
                
            except requests.exceptions.ConnectionError as e:
                error_msg = str(e)
                print(f"GraphQL request failed (attempt {attempt + 1}/{self.max_retries}): {error_msg}")
                
                if "Connection reset by peer" in error_msg or "Connection aborted" in error_msg:
                    print("  → Connection reset error detected")
                elif "timeout" in error_msg.lower():
                    print("  → Timeout error detected")
                elif "SSL" in error_msg or "certificate" in error_msg.lower():
                    print("  → SSL/certificate error detected")
                else:
                    print(f"  → Other connection error: {error_msg}")
                
                if attempt < self.max_retries - 1:
                    wait_time = self.exponential_backoff(attempt)
                    print(f"  → Waiting {wait_time:.1f} seconds before retry...")
                    time.sleep(wait_time)
                else:
                    print(f"  → Failed after {self.max_retries} attempts")
                    
            except requests.exceptions.Timeout as e:
                print(f"GraphQL request timeout (attempt {attempt + 1}/{self.max_retries}): {e}")
                if attempt < self.max_retries - 1:
                    wait_time = self.exponential_backoff(attempt)
                    print(f"  → Waiting {wait_time:.1f} seconds before retry...")
                    time.sleep(wait_time)
                else:
                    print(f"  → Failed after {self.max_retries} attempts")
                    
            except requests.exceptions.HTTPError as e:
                print(f"GraphQL HTTP error (attempt {attempt + 1}/{self.max_retries}): {e}")
                if e.response.status_code in [429, 500, 502, 503, 504]:
                    if attempt < self.max_retries - 1:
                        wait_time = self.exponential_backoff(attempt)
                        print(f"  → Waiting {wait_time:.1f} seconds before retry...")
                        time.sleep(wait_time)
                    else:
                        print(f"  → Failed after {self.max_retries} attempts")
                else:
                    print(f"  → Non-retryable HTTP error: {e.response.status_code}")
                    return None
                    
            except Exception as e:
                print(f"GraphQL request failed (attempt {attempt + 1}/{self.max_retries}): {e}")
                if attempt < self.max_retries - 1:
                    wait_time = self.exponential_backoff(attempt)
                    print(f"  → Waiting {wait_time:.1f} seconds before retry...")
                    time.sleep(wait_time)
                else:
                    print(f"  → Failed after {self.max_retries} attempts")
        
        return None
    
    def uniprot_from_pdb_chains_robust(self, pdb_chains: Set[str], batch_size: int = 50) -> pd.DataFrame:
        """Get UniProt IDs from PDB chains with robust error handling and batching"""
        
        pdb_chains_list = list(pdb_chains)
        total_chains = len(pdb_chains_list)
        print(f"Processing {total_chains} PDB chains in batches of {batch_size}")
        
        all_results = []
        successful_batches = 0
        failed_batches = 0
        
        # Process in batches
        for i in tqdm(range(0, total_chains, batch_size), desc="Processing PDB chain batches"):
            batch = pdb_chains_list[i:i + batch_size]
            batch_num = i // batch_size + 1
            
            print(f"\nProcessing batch {batch_num} ({len(batch)} chains)")
            
            # Create GraphQL query for this batch
            pdb_chains_string = '", "'.join(batch)
            query = f'''query {{
  polymer_entity_instances(instance_ids: ["{pdb_chains_string}"]) {{
    rcsb_id
    polymer_entity {{
      uniprots {{
        rcsb_id
      }}
    }}
  }}
}}'''
            
            # Make the request
            response_data = self.make_graphql_request(query)
            
            if response_data is not None:
                try:
                    # Process the response
                    batch_results = []
                    for record in response_data['data']['polymer_entity_instances']:
                        rcsb_id, chain = record['rcsb_id'].split('.')
                        if record['polymer_entity']['uniprots'] is not None:
                            uniprot_id = record['polymer_entity']['uniprots'][0]['rcsb_id']
                        else:
                            uniprot_id = ''
                        batch_results.append([rcsb_id, chain, uniprot_id])
                    
                    all_results.extend(batch_results)
                    successful_batches += 1
                    print(f"  ✓ Successfully processed {len(batch_results)} chains")
                    
                except Exception as e:
                    print(f"  → Error processing batch results: {e}")
                    failed_batches += 1
            else:
                print(f"  → Failed to retrieve data for batch {batch_num}")
                failed_batches += 1
            
            # Add a small delay between batches to be respectful to the server
            if i + batch_size < total_chains:
                time.sleep(1)
        
        print(f"\nSummary:")
        print(f"  Successful batches: {successful_batches}")
        print(f"  Failed batches: {failed_batches}")
        print(f"  Total results: {len(all_results)}")
        
        if all_results:
            return pd.DataFrame(all_results, columns=['pdb_id', 'chain', 'uniprot_id'])
        else:
            print("No results retrieved, returning empty DataFrame")
            return pd.DataFrame(columns=['pdb_id', 'chain', 'uniprot_id'])

# Convenience function
def uniprot_from_pdb_chains_robust(pdb_chains, batch_size=50):
    """Robust version of uniprot_from_pdb_chains with error handling and batching"""
    api = RobustPDBGraphQLAPI()
    return api.uniprot_from_pdb_chains_robust(pdb_chains, batch_size)

if __name__ == "__main__":
    # Test the robust function
    test_chains = {'1ABC.A', '2DEF.B', '3GHI.C'}
    print(f"Testing robust PDB GraphQL API with {len(test_chains)} chains")
    result = uniprot_from_pdb_chains_robust(test_chains, batch_size=10)
    print(f"Successfully retrieved data for {len(result)} chains")
    print(result)
