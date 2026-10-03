# Taguchi DOE Generator

A Streamlit app that builds randomized, strength-2 Taguchi orthogonal-array run sheets.

## Files

- `app.py` : Streamlit interface
- `taguchi_doe.py` : design engine (orthogonal arrays, DOF report, validation)
- `requirements.txt` : dependencies

## Run locally

```bash
pip install -r requirements.txt
streamlit run app.py
```

## Deploy on Streamlit Community Cloud

1. Push these files to a GitHub repository.
2. Go to https://share.streamlit.io and sign in with GitHub.
3. Choose **New app**, select the repository and branch, and set the main file to `app.py`.
4. Click **Deploy**.

## Notes

- Taguchi mode supports factors with 2, 3, 4, 5, 7, 8, 9, 11, 16... levels (prime powers). Factors can mix level counts.
- For 6, 10 or 12 levels, switch to full factorial.
- The design assigns main effects only. Leftover degrees of freedom are not automatically pure error.
