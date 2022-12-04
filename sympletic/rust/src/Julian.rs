mod Julian {
    struct Julian {
        date: f64
    }

    impl Julian {
        fn new(utc_datetime: DateTime<chrono::Utc>) -> Self {
            return Self { date: utc_datetime.timestamp()};
        }
        fn new() -> Self {
            return Self { date: 0};
        }
        
        fn get_date(&self) -> f64 {
            return self.date;
        }
    }
}