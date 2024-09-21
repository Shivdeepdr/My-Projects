SELECT * FROM actor;
SELECT min(payment_date) as min_payment, max(payment_date) as max_payment from payment;
SELECT * FROM payment;
SELECT * FROM payment ORDER BY payment_date DESC;

SELECT f.title, sum(p.amount) as revenue FROM payment as p
JOIN rental as r ON p.rental_id = r.rental_id
JOIN inventory as i ON r.inventory_id = i.inventory_id
JOIN film as f ON i.film_id = f.film_id
GROUP BY f.title
ORDER BY f.title ASC;

SELECT f.title, p.payment_id, p.customer_id, p.staff_id, p.amount, p.rental_id, r.inventory_id, r.rental_date, i.film_id FROM payment as p
JOIN rental as r ON p.rental_id = r.rental_id
JOIN inventory AS i ON r.inventory_id = i.inventory_id
JOIN film as f ON i.film_id = f.film_id;

SELECT f.title, ci.city,  sum(p.amount) as revenue FROM payment as p
JOIN rental as r ON p.rental_id = r.rental_id
JOIN inventory as i ON r.inventory_id = i.inventory_id
JOIN film as f ON i.film_id = f.film_id
JOIN customer as cu ON p.customer_id = cu.customer_id
JOIN address as a ON cu.address_id = a.address_id
JOIN city as ci ON a.city_id = ci.city_id
GROUP BY ci.city, f.title
ORDER BY revenue  DESC;

-- Creating star schema
-- Create dimdate
CREATE TABLE dimDate(
	dateKey INT NOT NULL PRIMARY KEY,
	date date NOT NULL,
	year SMALLINT NOT NULL,
	quarter SMALLINT NOT NULL,
	month SMALLINT NOT NULL,
	week SMALLINT NOT NULL,
	day SMALLINT NOT NULL,
	is_weekend BOOLEAN NOT NULL
);

SELECT * FROM information_schema.columns WHERE TABLE_NAME = 'dimdate'

INSERT INTO dimDate(datekey,date,year,quarter,month,week,day,is_weekend)
SELECT 
	distinct(to_char(payment_date ::DATE, 'yymmdd') :: integer ) as date_key,
	date(payment_date) as date,
	EXTRACT (year FROM payment_date) as year,
	EXTRACT (quarter FROM payment_date) as quarter,
	EXTRACT (month FROM payment_date) as month,
	extract (week from payment_date) as week,
	extract (day from payment_date) as day,
	case when extract(ISODOW from payment_date) in (6,7) THEN TRUE else FALSE end as is_weekend
FROM payment;

SELECT * FROM dimdate LIMIT(10);

-- Create dimcustomers

CREATE TABLE dimCustomers(
	customer_key SERIAL PRIMARY KEY,
	customer_id INT NOT NULL,
	first_name VARCHAR NOT NULL,
	last_name VARCHAR NOT NULL,
	email VARCHAR NOT NULL,
	address VARCHAR NOT NULL,
	address2 VARCHAR NOT NULL,
	district VARCHAR NOT NULL,
	city VARCHAR NOT NULL,
	country VARCHAR NOT NULL,
	postal_code VARCHAR NOT NULL,
	phone VARCHAR NOT NULL,
	active SMALLINT NOT NULL,
	create_date timestamp NOT NULL,
	start_date date NOT NULL,
	end_date date NOT NULL
);

SELECT * FROM information_schema.columns WHERE TABLE_NAME = 'dimcustomers'

INSERT INTO dimCustomers(customer_key,customer_id,first_name,last_name,email,
						  address,address2,district,city,country,postal_code,
						  phone,active,create_date,start_date,end_date)
						  
SELECT
c.customer_id as customer_key,
c.customer_id,
c.first_name,
c.last_name,
c.email,
a.address,
a.address2,
a.district,
ci.city,
co.country,
a.postal_code,
a.phone,
c.active,
c.create_date,
now() as start_date,
now() as end_date
FROM customer as c
JOIN address as a ON  c.address_id = a.address_id
JOIN city as ci ON ci.city_id = a.city_id
JOIN country as co ON ci.country_id = co.country_id;

SELECT * FROM dimcustomers LIMIT(15);

-- Create dimstore

CREATE TABLE dimStore(
	store_key SERIAL PRIMARY KEY,
	address VARCHAR NOT NULL,
	address2 VARCHAR,
	district VARCHAR NOT NULL,
	city VARCHAR NOT NULL,
	country VARCHAR NOT NULL,
	postal_code VARCHAR NOT NULL,
	managers_first_name VARCHAR NOT NULL,
	managers_last_name VARCHAR NOT NULL,
	start_date date NOT NULL,
	end_date date NOT NULL
);

SELECT * FROM information_schema.columns WHERE TABLE_NAME = 'dimstore' 

INSERT INTO dimStore(address,address2,district,city,country,postal_code,
					 managers_first_name,managers_last_name,start_date,end_date) 
SELECT
a.address as address,
a.address2 as address2,
a.district,
ci.city,
co.country,
a.postal_code,
s.first_name,
s.last_name,
now() as start_date,
now() as end_date
FROM address as a  
JOIN store as st ON a.address_id = st.address_id
JOIN city as ci ON a.city_id = ci.city_id
JOIN country as co ON ci.country_id = co.country_id
JOIN staff as s ON s.store_id = st.store_id;


SELECT * FROM dimstore;
SELECT address_id FROM address;

-- Create dimmovie

CREATE TABLE dimMovie(
	movie_key SERIAL PRIMARY KEY,
	film_id SMALLINT NOT NULL,
	title VARCHAR NOT NULL,
	description VARCHAR NOT NULL,
	release_year SMALLINT NOT NULL,
	length VARCHAR NOT NULL,
	ratings mpaa_rating NOT NULL,
	special_features VARCHAR NOT NULL
);

INSERT INTO dimMovie(movie_key,film_id,title,description,release_year,
					length,ratings,special_features)
SELECT
f.film_id as movie_key,
f.film_id,
f.title,
f.description,
f.release_year,
f.length,
f.rating,
f.special_features
FROM film as f;

SELECT * FROM dimmovie LIMIT(15);

-- create facttable

CREATE TABLE factTable(
	sales_key SERIAL PRIMARY KEY,
	date_key integer REFERENCES dimdate(datekey),
	customer_key integer REFERENCES dimCustomers(customer_key),
	movie_key integer REFERENCES dimMovie(movie_key),
	store_key integer REFERENCES dimstore(store_key),
	sales_amount numeric
);

SELECT * FROM information_schema.columns WHERE TABLE_NAME = 'facttable'

INSERT INTO  facttable(date_key,customer_key,movie_key,store_key,sales_amount)
SELECT 
	to_char(payment_date ::DATE, 'yymmdd') :: integer as date_key,
	p.customer_id as customer_key,
	dm.movie_key as movie_key,
	ds.store_key as store_key,
	p.amount as sales_amount
FROM payment as p
JOIN rental as r ON p.rental_id = r.rental_id
JOIN inventory as i ON r.inventory_id = i.inventory_id
JOIN dimmovie as dm ON i.film_id = dm.movie_key
JOIN staff as st ON r.staff_id = st.staff_id
JOIN dimstore as ds ON st.store_id = ds.store_key;


SELECT * FROM facttable;

SELECT ds.country, SUM(sales_amount) as revenue FROM facttable as ft
JOIN  dimmovie as dm ON (ft.movie_key = dm.movie_key)
JOIN dimstore as ds ON (ft.store_key = ds.store_key)
JOIN dimdate as dd ON (ft.date_key = dd.datekey)
GROUP BY (ds.country)
ORDER BY revenue ASC;


