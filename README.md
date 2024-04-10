# EntrezMetadataTools

ObtainEntrezMetadata.py3:

Capabilities:
- Obtain Entrez Metadata based on an arbitrary query. 
- Works irrespective of target Db (Bioproject, Biosample, Sra, etc).
- Built for high-throughput via batched queries. 
- Will pull IDs and metadata for 3M records in >2hrs with low overhead, **without parallelization.**
- Parallelization was added 2/27/2024 using `concurrent.futures`  `ThreadPoolExecutor` and `as_completed`. Script will now pull all Rna-Seq metadata from SRA in ~30 minutes
- XML parsing capability added 3/21/2024.
    - XML parser capable of serial flattening of compound nested XML structures.
    - XML parser capable of handling broken / incorrect XML formatted files of certain kinds.

2 tasks remain before versioning as 1.0:
1. Add Data Integration across multiple queries run using ObtainEntrezMetadata.py3. This will be with IntegrateEntrezMetadata.py3.
2. Add prioritization script to be able to deliver to Matt rank-ordered studies that are recommended.

Once the script provides helpful prioritization of studies, version as 1.0.

############## Notes pertaining to Data Merging 


Purpose of this Readme is to define terms used in the Metadata Cleaning and Merging Scripts (e.g. Metadata_Clean_and_Merge.py).

Formalized relationships used in relational database design are used to manage the metadata. We divide this into 2 procedures:

#####################################################################

Procedure 1: Pre-Merge Data-Processing of an individual table. This is, for example, what would be done to the raw data downloaded from NCBI (after
	processing from XML to XLSX, but before any additional cleaning. This is necessary because the raw data contains associated problems...

	Step 1.1: Frequently, columns will be present that are entirely or almost entirely empty. Here, we drop cols >95% empty (~85 of 138 for bioproject)
	Step 1.2: Copy the DF, then, on the copy, run lower() and strip() on every element to enhance identification of identity, bijectivity, etc
	Step 1.3: Identify any column that could conceivably serve as a primary key (strictly: is bijective unto itself) - (len (uniq(col)) == len(df) ).
	Step 1.4: Now, identify cols that are completely identical. remove from the duplicate and the original df.
	Step 1.5: Now, identify bijective cols. Remove one from the duplicate and original df, but create a foreign key
	Step 1.6: Identify surjective relationships. Use these to create Reference & Lookup Tables; STAR schemas (procedure further described below).

Procedure 2: Once 2 Processed DFs are obtained (e.g., one from bioproject and one from biosample), we wish to merge them. In so doing, we want to
apply many of the same principles again, though some (Step 1.1) arent needed and new steps need to be added. 
	Step 2.1: Identify any columns that could conceivably be used to join the two tables we heavily favor the use of identity for 2.1. Not abs requirement
	Step 2.2: Select a procedure for joining the two DBs. 
		If col1 == col2 			--> inner join.
		If col1 is a subset of col2 --> right join.
		If col2 is a subset of col1 -->  left join.
		Otherwise:					--> outer join.
	Step 2.3 - 2.5: After joining, repeat steps 1.4 - 1.6

#####################################################################
Terminology from Set Theory used in comments:

Injective Function (One-to-One Function): An injective function is one where every element of the domain maps to a unique element in the codomain. In other words, no two different elements in the domain map to the same element in the codomain. This ensures that the function has a kind of uniqueness for each input.

Surjective Function (Onto Function): A Surjective Function (crudely, one-to-many): is one where every element of the codomain is the image of at least one element from the domain. This means the function covers the entire codomain, leaving no element unmatched.

Bijective Function (One-to-One Correspondence): A bijective function combines the properties of both injective and surjective functions. It is both one-to-one and onto, meaning each element of the domain maps to a unique element in the codomain, and every element of the codomain is mapped by an element of the domain. Bijective functions are particularly important because they allow for inverse functions to be defined.

As examples of the above, consider:

# import pandas as pd
# import numpy as np

# Create a DataFrame with 8 unique values in column 'A'

    """ Example Usage """
    df = pd.DataFrame(
    {
    'A': np.arange(1, 9),               # Unique values from 1 to 8 (bijective)
    'B': np.random.permutation(np.arange(1, 9)),  # A permutation of A (bijective)
    'C': [1, 2, 2, 3, 3, 3, 4, 4],       # Repeated values, not unique (neither)
    'D': np.random.permutation(np.arange(1, 9)),  # Another permutation of A (bijective)
    'E': [1, 1, 2, 2, 3, 3, 4, 4],       # Subset of A's values but repeated (surjective (loosely depending on true domain def))
    'F': np.arange(1, 5).repeat(2)[:8],  # Unique values, partial coverage of A (injective)
    'G': np.arange(100, 108),            # Completely unique set (injective)
    'H': [1, 2, 3, 4, 5, 6, 7, 1]        # Duplicates but covers all A (surjective)
    } )
    print(df)

Note that we identify relationships of this kind between all pairs of columns in order to achieve greater data compression through the creation of ancillary tables.

#####################################################################################
Terminology from relational database design used in comments:

In relational database design, **primary keys**, **unique constraints**, and **foreign keys** are fundamental concepts used to ensure data integrity and establish relationships between tables. Here’s a detailed look at each of these components:

### Primary Keys
A **primary key** is a column (or a set of columns) in a database table that uniquely identifies each row in that table. The primary key must contain unique values, and it **** cannot contain null values **** (see below). This ensures that every record can be uniquely identified, which is crucial for relational integrity and efficient data retrieval.
- **Purpose**: To uniquely identify each record in a table.
- In this workflow, by limiting the search space to columns with 0% missingness and 100% uniqueness, we can rapidly identify candidates. NCBI Accession terms are preferred over any other kind of column if possible.

### Unique Constraints
A **unique constraint** ensures that all values in a column are different from each other. This constraint can apply to one or more columns, and unlike a primary key, it allows for null values (unless explicitly forbidden). A table can have multiple unique constraints.
- **Purpose**: To guarantee that no two rows have the same value in specific columns, thus maintaining uniqueness across the specified data.
- **Usage Example**: In a user table, you might have a unique constraint on the email address column to ensure that no two users can register with the same email address.

### Foreign Keys
A **foreign key** is a column (or collection of columns) in one table that uniquely identifies a row of another table or the same table (in case of a recursive relationship). The foreign key establishes a link between the data in two tables, typically enforcing a relationship where the foreign key data must match that in the primary key it references.
- **Purpose**: To create and enforce a link between the data in two tables, which is the cornerstone of the relational database model. This ensures the referential integrity of the data.
- **Usage Example**: In a database with tables for `Employees` and `Departments`, the `Employees` table could have a column `DepartmentID` that is a foreign key referencing the `Departments` table's primary key, ensuring that each employee is linked to a valid department.

### Lookup, Reference, and Star Table creation.
These elements help manage and maintain the consistency and accuracy of the data within relational databases. They are used to enforce data integrity, prevent invalid data entries, and help implement relationships between tables which are essential for complex queries and data analysis.

In this particular analysis, we identify surjective relationships, then remove col2, retaining col1, but we write only the unique relationships between col1 and col2.





