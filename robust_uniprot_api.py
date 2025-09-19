#!/usr/bin/env python3
"""
Robust UniProt API functions with comprehensive error handling
"""

import SPARQLWrapper
import itertools
import pandas as pd
import time
import random
from typing import List, Dict, Any
import sys

if sys.version_info >= (3, 12):
    from itertools import batched
else:
    try:
        from more_itertools import batched
    except ImportError:
        def batched(iterable, chunk_size):
            iterator = iter(iterable)
            while chunk := tuple(itertools.islice(iterator, chunk_size)):
                yield chunk

class RobustUniProtAPI:
    """Robust UniProt API client with retry logic and error handling"""
    
    def __init__(self, max_retries=5, base_delay=2, max_delay=60):
        self.max_retries = max_retries
        self.base_delay = base_delay
        self.max_delay = max_delay
        self.sparql_endpoint = "https://sparql.uniprot.org/sparql"
    
    def exponential_backoff(self, attempt: int) -> float:
        """Calculate exponential backoff delay with jitter"""
        delay = min(self.base_delay * (2 ** attempt), self.max_delay)
        jitter = random.uniform(0.1, 0.3) * delay
        return delay + jitter
    
    def make_sparql_request(self, query: str, timeout: int = 120) -> Dict[str, Any]:
        """Make a robust SPARQL request to UniProt"""
        
        for attempt in range(self.max_retries):
            try:
                sparql = SPARQLWrapper.SPARQLWrapper(self.sparql_endpoint)
                sparql.setQuery(query)
                sparql.setReturnFormat(SPARQLWrapper.JSON)
                sparql.setTimeout(timeout)
                
                # Add custom headers to mimic a real browser
                sparql.addCustomHttpHeader("User-Agent", 
                    "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36")
                sparql.addCustomHttpHeader("Accept", "application/sparql-results+json")
                sparql.addCustomHttpHeader("Accept-Language", "en-US,en;q=0.9")
                
                data = sparql.query().convert()
                return data
                
            except Exception as e:
                error_msg = str(e)
                print(f"SPARQL request failed (attempt {attempt + 1}/{self.max_retries}): {error_msg}")
                
                if "Connection reset by peer" in error_msg or "Connection aborted" in error_msg:
                    print("  → Connection reset error detected")
                elif "timeout" in error_msg.lower():
                    print("  → Timeout error detected")
                elif "urllib.error" in error_msg:
                    print("  → URL error detected")
                elif "allocate" in error_msg and "MB" in error_msg:
                    print("  → Memory allocation error detected - server out of memory")
                    print("  → This usually means the batch size is too large")
                    # For memory errors, don't retry immediately, return None
                    return None
                elif "QueryBadFormed" in error_msg or "badly formed" in error_msg.lower():
                    print("  → Query format error detected")
                    # For query format errors, don't retry
                    return None
                else:
                    print(f"  → Other error: {error_msg}")
                
                if attempt < self.max_retries - 1:
                    wait_time = self.exponential_backoff(attempt)
                    print(f"  → Waiting {wait_time:.1f} seconds before retry...")
                    time.sleep(wait_time)
                else:
                    print(f"  → Failed after {self.max_retries} attempts")
        
        return None
    
    def uniprot_protein_name(self, uniprot_ids: List[str], batch_size: int = 50) -> pd.DataFrame:
        """Get protein names from UniProt with robust error handling"""
        
        output = []
        total_ids = len(uniprot_ids)
        successful_batches = 0
        failed_batches = 0
        
        print(f"Processing {total_ids} UniProt IDs in batches of {batch_size}")
        
        for batch_num, uniprot_subset in enumerate(batched(uniprot_ids, batch_size), 1):
            uniprot_string = ' '.join([f'uniprotkb:{id}' for id in uniprot_subset])
            
            query_string = f"""
PREFIX up: <http://purl.uniprot.org/core/>
PREFIX uniprotkb: <http://purl.uniprot.org/uniprot/>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
SELECT DISTINCT ?uniprot_id ?reviewed ?name
WHERE
{{
  VALUES ?protein {{ {uniprot_string} }}
  BIND(substr(str(?protein), strlen(str(uniprotkb:))+1) AS ?uniprot_id)
  ?protein up:reviewed ?reviewed .
  OPTIONAL {{
    ?protein up:recommendedName ?recommended .
    ?recommended up:fullName ?name .
  }}
}}
"""
            
            print(f"Processing batch {batch_num} ({len(uniprot_subset)} IDs)")
            
            data = self.make_sparql_request(query_string)
            
            if data and "results" in data and "bindings" in data["results"]:
                batch_results = data["results"]["bindings"]
                print(f"  ✓ Successfully retrieved {len(batch_results)} results")
                
                for result in batch_results:
                    output.append({key: value['value'] for key, value in result.items()})
                
                successful_batches += 1
            else:
                print(f"  ✗ Failed to retrieve data for batch {batch_num}")
                failed_batches += 1
            
            # Add delay between batches
            if batch_num * batch_size < total_ids:
                time.sleep(2)
        
        print(f"\nSummary:")
        print(f"  Successful batches: {successful_batches}")
        print(f"  Failed batches: {failed_batches}")
        print(f"  Total results: {len(output)}")
        
        if output:
            return pd.DataFrame(output, columns=['uniprot_id', 'reviewed', 'name']).rename(
                columns={'uniprot_id': 'UniProt ID', 'reviewed': 'Reviewed', 'name': 'Protein Name'})
        else:
            return pd.DataFrame(columns=['UniProt ID', 'Reviewed', 'Protein Name'])
    
    def get_uniprot_site_annotations(self, uniprot_ids: List[str], batch_size: int = 50) -> pd.DataFrame:
        """Get site annotations from UniProt with robust error handling"""
        
        output = []
        total_ids = len(uniprot_ids)
        successful_batches = 0
        failed_batches = 0
        
        print(f"Processing {total_ids} UniProt IDs for site annotations in batches of {batch_size}")
        
        for batch_num, uniprot_subset in enumerate(batched(uniprot_ids, batch_size), 1):
            uniprot_string = ' '.join([f'uniprotkb:{id}' for id in uniprot_subset])
            
            query_string = f"""
PREFIX up: <http://purl.uniprot.org/core/>
PREFIX uniprotkb: <http://purl.uniprot.org/uniprot/>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX faldo: <http://biohackathon.org/resource/faldo#>
SELECT DISTINCT ?uniprot_id ?begin ?end ?site ?comment
WHERE
{{
    VALUES ?protein {{ {uniprot_string} }}
    BIND(substr(str(?protein), strlen(str(uniprotkb:))+1) AS ?uniprot_id)
  
    ?protein up:annotation ?annotation .
  {{ ?annotation a up:Binding_Site_Annotation }} UNION {{ ?annotation a up:Active_Site_Annotation }} .
    ?annotation rdf:type ?type .
    BIND(substr(str(?type), strlen(str(up:))+1) AS ?site)
    ?annotation up:range ?range .
    ?range faldo:begin/faldo:position ?begin .
    ?range faldo:end/faldo:position ?end .
    OPTIONAL
    {{
        ?annotation up:ligand ?ligand .
        ?ligand rdfs:comment ?comment .
    }}
}}
"""
            
            print(f"Processing batch {batch_num} ({len(uniprot_subset)} IDs)")
            
            data = self.make_sparql_request(query_string)
            
            if data and "results" in data and "bindings" in data["results"]:
                batch_results = data["results"]["bindings"]
                print(f"  ✓ Successfully retrieved {len(batch_results)} results")
                
                for result in batch_results:
                    output.append({key: value['value'] for key, value in result.items()})
                
                successful_batches += 1
            else:
                print(f"  ✗ Failed to retrieve data for batch {batch_num}")
                failed_batches += 1
            
            # Add delay between batches
            if batch_num * batch_size < total_ids:
                time.sleep(2)
        
        print(f"\nSummary:")
        print(f"  Successful batches: {successful_batches}")
        print(f"  Failed batches: {failed_batches}")
        print(f"  Total results: {len(output)}")
        
        if output:
            return pd.DataFrame(output, columns=['uniprot_id', 'site', 'begin', 'end', 'comment'])
        else:
            return pd.DataFrame(columns=['uniprot_id', 'site', 'begin', 'end', 'comment'])
    
    def get_uniprot_sequence(self, uniprot_ids: List[str], batch_size: int = 50) -> pd.DataFrame:
        """Get sequences from UniProt with robust error handling"""
        
        output = []
        total_ids = len(uniprot_ids)
        successful_batches = 0
        failed_batches = 0
        
        print(f"Processing {total_ids} UniProt IDs for sequences in batches of {batch_size}")
        
        for batch_num, uniprot_subset in enumerate(batched(uniprot_ids, batch_size), 1):
            uniprot_string = ' '.join([f'uniprotkb:{id}' for id in uniprot_subset])
            
            query_string = f"""
PREFIX up: <http://purl.uniprot.org/core/>
PREFIX uniprotkb: <http://purl.uniprot.org/uniprot/>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
SELECT DISTINCT ?uniprot_id ?sequence
WHERE
{{
    VALUES ?protein {{ {uniprot_string} }}
    BIND(substr(str(?protein), strlen(str(uniprotkb:))+1) AS ?uniprot_id)
    ?protein up:sequence ?isoform .
    ?isoform a up:Simple_Sequence ;
    rdf:value ?sequence .
}}
"""
            
            print(f"Processing batch {batch_num} ({len(uniprot_subset)} IDs)")
            
            data = self.make_sparql_request(query_string)
            
            if data and "results" in data and "bindings" in data["results"]:
                batch_results = data["results"]["bindings"]
                print(f"  ✓ Successfully retrieved {len(batch_results)} results")
                
                for result in batch_results:
                    output.append({key: value['value'] for key, value in result.items()})
                
                successful_batches += 1
            else:
                print(f"  ✗ Failed to retrieve data for batch {batch_num}")
                failed_batches += 1
            
            # Add delay between batches
            if batch_num * batch_size < total_ids:
                time.sleep(2)
        
        print(f"\nSummary:")
        print(f"  Successful batches: {successful_batches}")
        print(f"  Failed batches: {failed_batches}")
        print(f"  Total results: {len(output)}")
        
        if output:
            return pd.DataFrame(output, columns=['uniprot_id', 'sequence'])
        else:
            return pd.DataFrame(columns=['uniprot_id', 'sequence'])

# Convenience functions for backward compatibility
def uniprot_protein_name_robust(uniprot_ids, batch_size=100):
    """Robust version of uniprot_protein_name function"""
    api = RobustUniProtAPI()
    return api.uniprot_protein_name(uniprot_ids, batch_size)

def get_uniprot_site_annotations_robust(uniprot_ids, batch_size=100):
    """Robust version of get_uniprot_site_annotations function"""
    api = RobustUniProtAPI()
    return api.get_uniprot_site_annotations(uniprot_ids, batch_size)

def get_uniprot_sequence_robust(uniprot_ids, batch_size=100):
    """Robust version of get_uniprot_sequence function"""
    api = RobustUniProtAPI()
    return api.get_uniprot_sequence(uniprot_ids, batch_size)

if __name__ == "__main__":
    # Test the robust function
    test_ids = ['P12345', 'P67890', 'Q12345']
    print(f"Testing robust UniProt API with {len(test_ids)} IDs")
    result = uniprot_protein_name_robust(test_ids, batch_size=2)
    print(f"Successfully retrieved data for {len(result)} UniProt IDs")
