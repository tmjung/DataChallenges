#!/usr/bin/env python3
"""
Convert SQL dump to CSV format
Extracts data from phpMyAdmin SQL exports
"""

import re
import csv
import sys

def extract_columns_from_create_table(sql_content):
    """Extract column names from CREATE TABLE statement"""
    # Find the CREATE TABLE statement
    create_match = re.search(r'CREATE TABLE.*?\((.*?)\) ENGINE', sql_content, re.DOTALL)
    if not create_match:
        return None
    
    table_def = create_match.group(1)
    columns = []
    
    for line in table_def.split('\n'):
        line = line.strip()
        if line and not line.startswith('--') and not line.startswith('`id_'):
            # Extract column name (first word after backtick)
            match = re.match(r'`([^`]+)`', line)
            if match:
                columns.append(match.group(1))
    
    return columns

def parse_sql_values(values_str):
    """Parse a single VALUES tuple from SQL"""
    # Remove outer parentheses
    values_str = values_str.strip()[1:-1]
    
    values = []
    current = ""
    in_quotes = False
    escape_next = False
    
    for char in values_str:
        if escape_next:
            current += char
            escape_next = False
        elif char == '\\':
            escape_next = True
            current += char
        elif char == "'" and not escape_next:
            in_quotes = not in_quotes
            current += char
        elif char == ',' and not in_quotes:
            # This is a delimiter
            value = current.strip()
            # Remove surrounding quotes if present
            if value.startswith("'") and value.endswith("'"):
                value = value[1:-1]
                # Unescape quotes
                value = value.replace("\\'", "'").replace('\\\\', '\\')
            elif value == 'NULL':
                value = ''
            values.append(value)
            current = ""
        else:
            current += char
    
    # Don't forget the last value
    if current:
        value = current.strip()
        if value.startswith("'") and value.endswith("'"):
            value = value[1:-1]
            value = value.replace("\\'", "'").replace('\\\\', '\\')
        elif value == 'NULL':
            value = ''
        values.append(value)
    
    return values

def sql_dump_to_csv(sql_file, csv_file):
    """Convert SQL dump to CSV"""
    
    with open(sql_file, 'r', encoding='utf-8') as f:
        sql_content = f.read()
    
    # Extract column names
    columns = extract_columns_from_create_table(sql_content)
    if not columns:
        print("Error: Could not extract column names from CREATE TABLE")
        return False
    
    print(f"Found {len(columns)} columns: {columns[:5]}...")
    
    # Extract all INSERT statements
    insert_matches = re.finditer(r'INSERT INTO.*?\((.*?)\) VALUES\s*(.*?)(?=;|INSERT)', sql_content, re.DOTALL)
    
    all_rows = []
    for insert_match in insert_matches:
        # Extract the columns being inserted
        cols_str = insert_match.group(1)
        insert_cols = [c.strip().strip('`') for c in cols_str.split(',')]
        
        # Extract values
        values_str = insert_match.group(2).strip()
        
        # Split by ),(  to separate multiple rows
        row_matches = re.finditer(r'\(([^()]*(?:\([^()]*\)[^()]*)*)\)', values_str)
        
        for row_match in row_matches:
            try:
                values = parse_sql_values('(' + row_match.group(1) + ')')
                
                # Create row dict with all columns
                row_dict = {col: '' for col in columns}
                
                # Fill in the values for the inserted columns
                for col, val in zip(insert_cols, values):
                    if col in row_dict:
                        row_dict[col] = val
                
                all_rows.append(row_dict)
            except Exception as e:
                print(f"Warning: Could not parse row: {e}")
                continue
    
    print(f"Extracted {len(all_rows)} rows")
    
    # Write to CSV
    try:
        with open(csv_file, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=columns)
            writer.writeheader()
            writer.writerows(all_rows)
        print(f"Successfully wrote {len(all_rows)} rows to {csv_file}")
        return True
    except Exception as e:
        print(f"Error writing CSV: {e}")
        return False

if __name__ == "__main__":
    sql_file = "/home/users8/bis/s7800363/Dev/DataChallenges/sources/ffm_vfpa_eisenzeit"
    csv_file = "/home/users8/bis/s7800363/Dev/DataChallenges/data/ffm_vfpa_eisenzeit.csv"
    
    print(f"Converting {sql_file} to {csv_file}...")
    if sql_dump_to_csv(sql_file, csv_file):
        print("Conversion complete!")
    else:
        print("Conversion failed!")
        sys.exit(1)
