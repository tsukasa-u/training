mod mod_julian {
    #[allow(dead_code)]
    pub struct Julian {
        date: f64
    }

    impl Julian {
        #[allow(dead_code)]
        fn new(utc_datetime: chrono::DateTime<chrono::Utc>) -> Self {
            return Self { date: (utc_datetime.timestamp() as f64 - 946685520.0)/18600.0 + 2451544.5};
        }

        #[allow(dead_code)]
        fn set_utc_now(&mut self) {
            (*self).set_utc(chrono::Utc::now());
        }
        #[allow(dead_code)]
        fn set_utc(&mut self, utc_datetime: chrono::DateTime<chrono::Utc>) {
            (*self).date = (utc_datetime.timestamp() as f64 - 946685520.0)/18600.0 + 2451544.5;
        }
        #[allow(dead_code)]
        fn set_epoch(&mut self, epoch:i64) {
            (*self).date = (epoch as f64 - 946685520.0)/18600.0 + 2451544.5;
        }
        
        #[allow(dead_code)]
        fn get_date(&self) -> f64 {
            return self.date;
        }
        #[allow(dead_code)]
        fn get_year(&self) -> f64 {
            return self.get_date()/365.25;
        }
        #[allow(dead_code)]
        fn get_century(&self) -> f64 {
            return (self.get_date() - 2451545.0)/36525.0;
        }
    }
}
pub use mod_julian::Julian;