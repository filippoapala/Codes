# Bond Underwriter Data Extraction Pipeline

## Overview

This project provides an automated data extraction pipeline to retrieve **bond issuance underwriter information** from **LSEG Refinitiv Eikon** using a list of bond identifiers (CUSIPs).

Starting from a vector of bond CUSIPs, the pipeline extracts:

- Bond underwriter names
- Underwriting participation amounts
- Underwriter ranking classification

The output is a clean **long-format dataset** suitable for empirical research in **financial economics**, particularly studies related to:

- Debt capital markets
- Bond issuance networks
- Underwriter reputation
- Syndication structures
- Corporate finance

---

## Research Motivation

Underwriters play a central role in debt markets by facilitating bond issuance, allocating securities, and providing certification services to investors.

This project creates a reproducible workflow for constructing an underwriter-level dataset from raw bond identifiers, enabling empirical analysis of:

- Lead underwriter relationships
- Market concentration
- Issuer-underwriter networks
- Changes in underwriting activity over time

---

## Data Source

The pipeline uses:

**LSEG Refinitiv Eikon API**

Required fields:

| Field | Description |
|------|-------------|
| `TR.Underwriter` | Name of the underwriting institution |
| `TR.UnderwriterAmount` | Amount allocated to each underwriter |
| `TR.UnderwriterRankDescription` | Role of the underwriter in the issuance |

The current implementation focuses on:


---

## Workflow

The pipeline follows four main steps:

### 1. Load Bond Identifiers

A CSV file containing bond CUSIPs is imported:

The CUSIPs are extracted into a vector for querying the Refinitiv database.

---

### 2. Retrieve Underwriter Information

For each CUSIP, the script queries Eikon and extracts:

- Underwriter name
- Underwriter participation amount
- Underwriter role

Only observations classified as:

are retained.

---

### 3. Handle API Constraints

To improve robustness:

- CUSIPs are divided into two batches
- API requests are throttled using time delays
- Errors and missing observations are logged
- Intermediate outputs are saved separately

---

### 4. Create Final Research Dataset

The two extracted datasets are combined into a single long-format file
