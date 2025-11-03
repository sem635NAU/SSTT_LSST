import pandas as pd

def make_csv(dictionary_of_asteroids, obs_or_fut, start_date, end_date):
    da = dictionary_of_asteroids
    df = pd.DataFrame(da, index=da['ssnamenr'])
    df = df.sort_index()
    df.to_csv(f'LSST_{obs_or_fut}_{start_date}_{end_date}.csv', index=False)

def auto_email(csv_file):
    import smtplib
    from email.mime.text import MIMEText
    from email.mime.multipart import MIMEMultipart
    