import glob
import os

import pandas as pd
import xmltodict
import tarfile


def parse_allosteric_site(alloteric_site):
    """Convert the list of residues in ASD from "Chain A:HIS25,TYR258; Chain B:VAL325" to ['A-HIS-25', 'A-TYR-258', 'B-VAL-325']"""
    residues = []
    for chain_string in alloteric_site.split('; '):
        chain_name, residue_string = chain_string.split(':')
        chain_id = chain_name[-1]
        for residue in residue_string.split(','):
            res_name, res_id = residue[:3], residue[3:]
            residues.append(f'{chain_id}-{res_name}-{res_id}')
    return residues


def parse_asd_xml(xml_string: str) -> tuple:
    """parse a directory containing XML files from the ASD database to create a pandas DataFrame"""
    xml_string = xml_string.replace('&#x2;', '')  # Remove invalid XML character
    protein = xmltodict.parse(xml_string)

    target_gene = ''
    organism_record = protein['Organism_Record']
    target_id = organism_record['Organism_ID']

    if 'Gene' in organism_record and 'Gene_Name' in organism_record['Gene']:
        target_gene = organism_record['Gene']['Gene_Name']

    organism = organism_record['Organism']
    molecule_class = organism_record['Molecule_Class']

    # Get EC Number
    if 'Enzyme_Nomenclature_List' in organism_record:
        enzyme_nomenclature = organism_record['Enzyme_Nomenclature_List']['Enzyme_Nomenclature']
        if isinstance(enzyme_nomenclature, list):
            ec_number = []
            for enzyme in enzyme_nomenclature:
                ec_number.append(enzyme['Enzyme_DB_ID'])
        else:
            ec_number = [enzyme_nomenclature['Enzyme_DB_ID']]
    else:
        ec_number = []
    
    # Check if the list of allosteric sites in present in the XML file
    if 'Allosteric_Site_List' in organism_record:
        allosteric_sites = organism_record['Allosteric_Site_List']['Allosteric_Site']

        # Check if more than one allosteric site is present
        if isinstance(allosteric_sites, list):
            allosteric_sites = allosteric_sites
        else:
            allosteric_sites = [allosteric_sites]
    else:
        allosteric_sites = []

    output = []
    for site in allosteric_sites:
        if site:
            pdb_uniprot = site['PDB_UniProt_ID'].upper() if 'PDB_UniProt_ID' in site else ''
            allosteric_pdb = site['Allosteric_PDB'].upper()
            modulator_serial = site['Modulator_ASD_ID']
            modulator_alias = site['Modulator_Alias']
            modulator_chain = site['Modulator_Chain']
            modulator_class = site['Modulator_Class'] if 'Modulator_Class' in site else ''
            modulator_feature = site['Modulator_Feature'] if 'Modulator_Feature' in site else ''
            modulator_name = site['Modulator_Name'] if 'Modulator_Name' in site else ''
            modulator_resi = site['Modulator_Residue'] if 'Modulator_Residue' in site else ''
            function = site['Function'] if 'Function' in site else ''
            position = site['Position'] if 'Position' in site else ''
            pubmed_id = site['PubMed_ID'] if 'PubMed_ID' in site else ''
            ref_title = site['PubMed_Title'] if 'PubMed_Title' in site else ''
            site_overlap = site['Site_Overlap'] if 'Site_Overlap' in site else ''
            allosteric_site_residue = parse_allosteric_site(site['Allosteric_Site_Residue']) if 'Allosteric_Site_Residue' in site else []

            output.append([
                target_id, target_gene, organism, pdb_uniprot, allosteric_pdb, molecule_class, ec_number,
                modulator_serial, modulator_alias, modulator_chain,
                modulator_class, modulator_feature, modulator_name,
                modulator_resi, function, position, pubmed_id, ref_title,
                site_overlap, allosteric_site_residue])
    return output


def asd_to_df(asd_xml_archive):
    asd_data = []
    with tarfile.open(asd_xml_archive, 'r:gz') as tar_data:
        for xml_file in tar_data:
            xml_data = tar_data.extractfile(xml_file)
            if xml_data is not None:
                xml_string = xml_data.read().decode("utf-8")
                asd_data.extend(parse_asd_xml(xml_string))

    return pd.DataFrame(asd_data, columns=[
        'Protein ASD ID',
        'Gene',
        'Organism',
        'UniProt ID',
        'PDB ID',
        'Protein Class',
        'EC Number',
        'Modulator ASD ID',
        'Modulator Alias',
        'Modulator Chain',
        'Modulator Class',
        'Allosteric Activity',
        'Modulator Name',
        'Modulator Residue ID',
        'ASD Function',
        'Position',
        'PubMed',
        'Reference Title',
        'Site Overlap',
        'ASD Allosteric Site Residues'])


if __name__ == "__main__":
    
    """Download ASD_Release_202306_XF.tar.gz from
    https://mdl.shsmu.edu.cn/ASD/module/download/download.jsp?tabIndex=1"""
    
    df_asd = parse_asd_xml('../ASD_Release_202306_XF.tar.gz')
    df_asd.to_csv('../ASD_Release_202306.csv', index=False)
