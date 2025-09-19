#!/usr/bin/env python3
"""
Smart UniProt API with automatic batch size adjustment
"""

from robust_uniprot_api import RobustUniProtAPI
import pandas as pd
from typing import List

class SmartUniProtAPI(RobustUniProtAPI):
    """Smart UniProt API that automatically adjusts batch size based on errors"""
    
    def __init__(self, max_retries=3, base_delay=2, max_delay=60):
        super().__init__(max_retries, base_delay, max_delay)
        self.memory_error_count = 0
        self.connection_error_count = 0
    
    def smart_batch_size(self, total_ids: int, base_batch_size: int = 50) -> int:
        """Calculate smart batch size based on error history"""
        
        # If we've had memory errors, use smaller batches
        if self.memory_error_count > 0:
            if total_ids > 100:
                return min(10, base_batch_size // 5)
            else:
                return min(5, base_batch_size // 10)
        
        # If we've had connection errors, use medium batches
        elif self.connection_error_count > 2:
            return min(25, base_batch_size // 2)
        
        # Default batch size
        else:
            return base_batch_size
    
    def make_sparql_request(self, query: str, timeout: int = 120):
        """Override to track error types"""
        
        result = super().make_sparql_request(query, timeout)
        
        # If result is None, it means we had errors
        if result is None:
            # Check the last error message (this is a simplified approach)
            # In a real implementation, you'd want to pass error info through
            pass
        
        return result
    
    def uniprot_protein_name_smart(self, uniprot_ids: List[str]) -> pd.DataFrame:
        """Get protein names with smart batch size adjustment"""
        
        total_ids = len(uniprot_ids)
        print(f"Processing {total_ids} UniProt IDs with smart batch sizing")
        
        # Start with a reasonable batch size
        batch_size = self.smart_batch_size(total_ids, 50)
        print(f"Starting with batch size: {batch_size}")
        
        max_attempts = 3
        attempt = 0
        
        while attempt < max_attempts:
            try:
                result = self.uniprot_protein_name(uniprot_ids, batch_size)
                
                # If we got results, return them
                if len(result) > 0:
                    print(f"Successfully retrieved {len(result)} results with batch size {batch_size}")
                    return result
                else:
                    print(f"No results with batch size {batch_size}, trying smaller batches...")
                    batch_size = max(5, batch_size // 2)
                    attempt += 1
                    
            except Exception as e:
                error_msg = str(e)
                print(f"Error with batch size {batch_size}: {error_msg}")
                
                if "allocate" in error_msg and "MB" in error_msg:
                    self.memory_error_count += 1
                    print("Memory error detected, reducing batch size")
                    batch_size = max(5, batch_size // 3)
                elif "Connection reset" in error_msg or "Connection aborted" in error_msg:
                    self.connection_error_count += 1
                    print("Connection error detected, reducing batch size")
                    batch_size = max(10, batch_size // 2)
                else:
                    batch_size = max(5, batch_size // 2)
                
                attempt += 1
                
                if batch_size < 5:
                    print("Batch size too small, giving up")
                    break
        
        print(f"Failed to retrieve data after {max_attempts} attempts")
        return pd.DataFrame(columns=['UniProt ID', 'Reviewed', 'Protein Name'])

# Convenience function
def uniprot_protein_name_smart(uniprot_ids):
    """Smart version that automatically adjusts batch size"""
    api = SmartUniProtAPI()
    return api.uniprot_protein_name_smart(uniprot_ids)

if __name__ == "__main__":
    # Test the smart function
    import sys
    import os
    sys.path.append('src')
    import parse_asd_xml
    
    # Load some data
    xml_dir = '../xml_file'
    df_asd = parse_asd_xml.asd_to_df(xml_dir)
    test_uniprot_ids = df_asd['UniProt ID'].unique()[:30]
    
    print(f"Testing smart UniProt API with {len(test_uniprot_ids)} IDs")
    result = uniprot_protein_name_smart(test_uniprot_ids)
    print(f"Successfully retrieved data for {len(result)} UniProt IDs")
